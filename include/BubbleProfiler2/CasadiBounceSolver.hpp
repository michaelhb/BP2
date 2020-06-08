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
    double V0;
    std::vector<double> Phi0;
    std::vector<double> U0;
};

struct NLP {
    casadi::Function nlp;
    
    // Separate T/V for ansatz / return is 
    // ugly and should be done better 
    casadi::Function T_a;
    casadi::Function T_ret;
    casadi::Function V_a;
    casadi::Function V_ret;
    casadi::Function Phi_ret;
};

class CasadiBounceSolver {

public:
    CasadiBounceSolver(casadi::Function potential_, int n_phi_,
        casadi::SXVector v_params_ = {}, int n_dims_ = 3, int N_ = 100);

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
    double grid_scale = 15.0; // Multiplying factor for gamma (TODO make dynamic)
    NLP nlp; // Algebraic representation of optimisation problem
    std::vector<std::vector<double> > C; // Coefficients of the collocation equation
    std::vector<double> D; // Coefficients of the continuity equation
    std::vector<double> B; // Coefficients for Gaussian quadrature 
    std::vector<casadi::Polynomial> P; // Basis of legendre polynomials

    //! Construct the parametric NLP we solve to find the bounce
    NLP get_nlp(const casadi::Function& potential) const;

    //! Time at element k, collocation point j
    double t_kj(int k, int j) const {
        return t_k[k] + h_k[k]*tau_root[j];
    }

    //! Transformation from compact to semi infinite domain
    double Tr(double tau) const {
        return grid_scale*log(2.0/(1.0 - tau));
    }

    //! Derivative of monotonic transformation
    double Tr_dot(double tau) const {
        return grid_scale / (1.0 - tau);
    }

    //! Ansatz in semi-infinite coordinates
    casadi::DM ansatz(double rho, casadi::DM true_vac, casadi::DM false_vac, 
        double r0, double sigma) const {        
        return true_vac + 0.5*(false_vac - true_vac)*(1 
            + tanh((rho - r0) / sigma)
            + exp(-rho)/(sigma*std::pow(cosh(r0/sigma),2)));
    }

    //! Derivative of ansatz in semi-infinite coordinates
    casadi::DM ansatz_dot(double rho, casadi::DM true_vac, casadi::DM false_vac,
        double r0, double sigma) const {
        return ((false_vac - true_vac)/2.0)*(
                1.0/(sigma*std::pow(cosh((rho-r0)/sigma), 2)) -
                exp(-rho)/(sigma*std::pow(cosh(r0/sigma), 2)));
    }

    //! Find ansatz parameters (r0, sigma) such that V[phi] = V0_target
    Ansatz get_ansatz(casadi::Function fV, 
            std::vector<double> grid_pars, 
            std::map<std::string, double> v_pars, 
            casadi::DM true_vac, casadi::DM false_vac) const;

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

    //! Concatenate vector of parametric input values
    std::vector<double> get_grid_pars() const {
        // This is a standin to allow for more flexibility later
        std::vector<double> par0;
        par0.insert(par0.end(), h_k.begin(), h_k.end());

        for (int k = 0; k < N; ++k) {
            for (int j = 0; j <= d; ++j) {
                par0.push_back(Tr(t_kj(k,j)));
            }
        }

        for (int k = 0; k < N; ++k) {
            for (int j = 0; j <= d; ++j) {
                par0.push_back(Tr_dot(t_kj(k,j)));
            }
        }
        return par0;
    }

    
};

}

#endif
