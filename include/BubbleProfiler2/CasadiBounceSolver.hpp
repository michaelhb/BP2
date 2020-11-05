#ifndef BUBBLEPROFILER2_CASADIBOUNCESOLVER_HPP_INCLUDED
#define BUBBLEPROFILER2_CASADIBOUNCESOLVER_HPP_INCLUDED

#include <BubbleProfiler2/BouncePath.hpp>

#include <memory>
#include <chrono>
#include <exception>
#include <math.h> 
#include <casadi/casadi.hpp>

namespace BubbleProfiler2 {

struct Ansatz {
    // TODO: omg make this a class already

    double V0; // (actual) value of V constraint
    double r0; // Estimated wall location (on [0,inf])
    double sigma; // Estimate wall thickness scale
    casadi::DM true_vac;
    casadi::DM false_vac;
    std::vector<double> Phi0; // Ansatz field values
    std::vector<double> U0; // Anstatz control values
};

struct NLP {
    casadi::Function nlp;
    casadi::Function ansatz_nlp;

    // Separate T/V for ansatz / return is 
    // ugly and should be done better 
    casadi::Function T_a; // T function (for ansatz)
    casadi::Function T_ret; // T function (for return)
    casadi::Function V_a; // V function (for ansatz)
    casadi::Function V_ret; // T function (for return)
    casadi::Function Phi_ret; // Extract state variables from solution
};

struct Static_bounds {
    // Control variable bounds
    std::vector<double> lbU;
    std::vector<double> ubU;

    // State variable bounds
    std::vector<double> lbPhi; 
    std::vector<double> ubPhi;

    // Constraint function bounds
    std::vector<double> lbg;
    std::vector<double> ubg;
};

//! Data class to represent information about the 
// compactified spatial grid.
class CompactGrid {
public:
    CompactGrid(
        double scale_,
        std::vector<double> h_k_,
        std::vector<double> gamma_kj_,
        std::vector<double> gammadot_kj_) : 
        h_k(h_k_), gamma_kj(gamma_kj_), gammadot_kj(gammadot_kj_), scale(scale_){

        concat = std::vector<double>(h_k);
        concat.insert(concat.end(), gamma_kj.begin(), gamma_kj.end());
        concat.insert(concat.end(), gammadot_kj.begin(), gammadot_kj.end());        
    }

    //! Scale factor for this grid
    double grid_scale() {
        return scale;
    }

    //! Finite element widths
    std::vector<double> get_h_k() {
        return h_k;
    }

    //! Monotonic transformation evaluated on all grid points
    std::vector<double> get_gamma_kj() {
        return gamma_kj;
    }

    //! Derivative of monotonic transformation evaluated on all grid points
    std::vector<double> get_gammadot_kj() {
        return gammadot_kj;
    }

    //! Concatenate in canonical order for NLP input
    std::vector<double> concatenate() {
        return concat;
    }

private:
    double scale;
    std::vector<double> h_k;
    std::vector<double> gamma_kj;
    std::vector<double> gammadot_kj;
    std::vector<double> concat;
};

class CasadiBounceSolver {

public:
    CasadiBounceSolver(casadi::Function potential_, int n_phi_,
        casadi::SXVector v_params_ = {}, int n_dims_ = 3, int N_ = 100, bool quiet_ = false);

    BouncePath solve(const std::vector<double>& true_vacuum, const std::vector<double>& false_vacuum,
        std::map<std::string, double> v_pars) const;

private:
    casadi::Function potential;
    casadi::SXVector v_params; // Vector of auxiliary parameter symbols
    int n_phi; // Number of field dimensions
    int n_dims; // Number of spatial dimensions
    double S_n; // Surface area of (d-1)-sphere
    int N; // Number of finite elements
    int d = 3; // Degree of interpolating polynomials
    std::vector<double> tau_root; // Vector of collocation points on [0,1)
    std::vector<double> t_k; // Element start times
    std::vector<double> h_k; // Element widths
    double default_grid_scale = 15.0; // Multiplying factor for gamma (TODO make dynamic)
    double r_point_frac = 0.7; // Fraction of points to place before the bubble wall when solving
    NLP nlp; // Algebraic representation of optimisation problem
    std::vector<std::vector<double> > C; // Coefficients of the collocation equation
    std::vector<double> D; // Coefficients of the continuity equation
    std::vector<double> B; // Coefficients for Gaussian quadrature 
    std::vector<casadi::Polynomial> P; // Basis of legendre polynomials
    bool quiet; // Suppress IPOPT output

    //! Construct the parametric NLP we solve to find the bounce
    NLP get_nlp(const casadi::Function& potential) const;

    //! Get static bounds on states, controls and constraints
    Static_bounds get_static_bounds(std::vector<double> false_vacuum) const;

    //! Time at element k, collocation point j
    double t_kj(int k, int j) const {
        return t_k[k] + h_k[k]*tau_root[j];
    }

    //! Transformation from compact to semi infinite domain
    template <typename T>
    T Tr(double tau, T grid_scale) const {
        using namespace casadi;
        return grid_scale*log(2.0/(1.0 - tau));
    }

    //! Derivative of monotonic transformation
    template <typename T>
    T Tr_dot(double tau, T grid_scale) const {
        return grid_scale / (1.0 - tau);
    } 

    //! Ansatz in semi-infinite coordinates
    casadi::DM ansatz(double rho, casadi::DM true_vac, casadi::DM false_vac, 
        double r0, double sigma) const {        
        return true_vac + 0.5*(false_vac - true_vac)*(1 
            + tanh((rho - r0) / sigma)
            + exp(-rho)/(sigma*std::pow(cosh(r0/sigma),2)));
    }

    //! Symbolic expansion of ansatz at Tr(t_kj)
    casadi::SX ansatz(int k, int j, casadi::SX true_vac, casadi::SX false_vac, 
        casadi::SX r0, casadi::SX sigma, casadi::SXVector gamma_par) const {

        using namespace casadi;
        SX rho = gamma_par[k](j);

        return true_vac + 0.5*(false_vac - true_vac)*(1 
            + tanh((rho - r0) / sigma)
            + (exp(-rho))/((sigma*pow(cosh(r0/sigma),2))));
    }
        
    //! Derivative of ansatz in semi-infinite coordinates
    casadi::DM ansatz_dot(double rho, casadi::DM true_vac, casadi::DM false_vac,
        double r0, double sigma) const {
        return ((false_vac - true_vac)/2.0)*(
                1.0/(sigma*std::pow(cosh((rho-r0)/sigma), 2)) -
                exp(-rho)/(sigma*std::pow(cosh(r0/sigma), 2)));
    }

    //! Symbolic expansion of ansatz derivative at Tr(t_kj)
    casadi::SX ansatz_dot(int k, int j,  casadi::SX true_vac, casadi::SX false_vac, 
        casadi::SX r0, casadi::SX sigma, casadi::SXVector gamma_par) const {

        using namespace casadi;
        SX rho = gamma_par[k](j);

        return ((false_vac - true_vac)/2.0)*(
                1.0/(sigma*pow(cosh((rho-r0)/sigma), 2)) -
                exp(-rho)/(sigma*pow(cosh(r0/sigma), 2)));
    }

    //! Find ansatz parameters (r0, sigma) such that V[phi] = V0_target
    Ansatz get_initial_ansatz(casadi::Function fV, casadi::Function fT,
            CompactGrid& grid,
            std::map<std::string, double> v_pars, 
            casadi::DM true_vac, casadi::DM false_vac) const;

    //! Find the grid scale that places r_point_frac points before 
    // the estimated bubble wall
    double find_grid_scale(double r0) const {
        double t_point_frac = 2.0*r_point_frac - 1;
        return r0/(log((2.0)/(1 - t_point_frac)));
    }

    //! Rescale the ansatz solution, after adjusting the grid scale
    Ansatz rescale_ansatz(Ansatz ansatz, CompactGrid& grid) const;

    //! Given r0, sigma, and a CompactGrid, get the corresponding Ansatz
    // TODO: Make Ansatz a class already, this is terrible
    Ansatz get_ansatz(double r0, double sigma, double V0, std::vector<double> true_vac, 
        std::vector<double> false_vac, CompactGrid grid) const;

    //! Utility method: generate subscripted variable names
    std::string varname(std::string prefix, std::vector<int> indices) const {
        std::ostringstream ss;
        ss << prefix;
        for (auto &ix : indices) {
            ss << "_" << ix;
        }
        return ss.str();
    }

    //! Utility method: populate SXDict with parameter arguments
    void add_param_args(casadi::SXDict& dict) const {
        for (int i = 0; i < v_params.size(); ++i) {
            dict[v_params[i].name()] = v_params[i];
        }
    }

    //! Utility method: append one vector to another
    void append_d(std::vector<double> &to, std::vector<double> from) const {
        // Weird function signature is because I want to pass rValues for the 
        // second argument. There's probably a less stupid way to get what I 
        // want...
        to.insert(to.end(), from.begin(), from.end());
    }

    //! Get grid information for a given grid scale
    CompactGrid get_grid(double grid_scale) const {
        std::vector<double> gamma_kj, gammadot_kj;

        for (int k = 0; k < N; ++k) {
            for (int j = 0; j <= d; ++j) {
                gamma_kj.push_back(Tr(t_kj(k, j), grid_scale));
                gammadot_kj.push_back(Tr_dot(t_kj(k, j), grid_scale));
            }
        }

        return CompactGrid(grid_scale, h_k, gamma_kj, gammadot_kj);
    }
};

}

#endif
