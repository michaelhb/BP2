#ifndef BUBBLEPROFILER2_POLYNOMIALS_HPP_INCLUDED
#define BUBBLEPROFILER2_POLYNOMIALS_HPP_INCLUDED

#include <BubbleProfiler2/CasadiBounceSolver.hpp>
#include <BubbleProfiler2/BouncePath.hpp>
#include <BubbleProfiler2/CasadiFindMinimum.hpp>
#include <cassert>

#define sqr(x) (x)*(x)

struct PolynomialPotential {
    casadi::Function fV;
    casadi::SX phi;
    casadi::SX delta;
};

std::vector<std::vector<double>> coeffs = {
    {1.8, 0.2},
    {0.684373, 0.181928, 0.295089},
    {0.534808, 0.77023, 0.838912, 0.00517238},
    {0.4747, 0.234808, 0.57023, 0.138912, 0.517238},
    {0.34234, 0.4747, 0.234808, 0.57023, 0.138912, 0.517238},
    {0.5233, 0.34234, 0.4747, 0.234808, 0.57023, 0.138912, 0.517238},
    {0.2434, 0.5233, 0.34234, 0.4747, 0.234808, 0.57023, 0.138912, 0.51723},
    {0.21, 0.24, 0.52, 0.34, 0.47, 0.23, 0.57, 0.14, 0.52},
    {0.12, 0.21, 0.24, 0.52, 0.34, 0.47, 0.23, 0.57, 0.14, 0.52},
    {0.23, 0.21, 0.21, 0.24, 0.52, 0.34, 0.47, 0.23, 0.57, 0.14, 0.52},
    {0.12, 0.11, 0.12, 0.21, 0.24, 0.52, 0.34, 0.47, 0.23, 0.57, 0.14, 0.52},
    {0.54, 0.47, 0.53, 0.28, 0.35, 0.27, 0.42, 0.59, 0.33, 0.16, 0.38, 0.35, 0.17},
    {0.39, 0.23, 0.26, 0.40, 0.11, 0.42, 0.41, 0.27, 0.42, 0.54, 0.18, 0.59, 0.13, 0.29},
    {0.21, 0.22, 0.22, 0.23, 0.39, 0.55, 0.43, 0.12, 0.16, 0.58, 0.25, 0.50, 0.45, 0.35, 0.45},
    {0.42, 0.34, 0.43, 0.22, 0.59, 0.41, 0.58, 0.41, 0.26, 0.45, 0.16, 0.31, 0.39, 0.57, 0.43, 0.10},
    {0.24, 0.35, 0.39, 0.56, 0.37, 0.41, 0.52, 0.31, 0.52, 0.22, 0.58, 0.39, 0.39, 0.17, 0.46, 0.30, 0.37},
    {0.18, 0.17, 0.30, 0.22, 0.38, 0.48, 0.11, 0.49, 0.43, 0.47, 0.21, 0.29, 0.32, 0.36, 0.30, 0.56, 0.46, 0.42},
    {0.40, 0.14, 0.10, 0.43, 0.39, 0.27, 0.33, 0.59, 0.48, 0.36, 0.24, 0.28, 0.51, 0.59, 0.40, 0.39, 0.24, 0.35, 0.20},
    {0.42, 0.11, 0.47, 0.13, 0.16, 0.24, 0.58, 0.53, 0.38, 0.44, 0.18, 0.46, 0.47, 0.27, 0.53, 0.24, 0.33, 0.40, 0.32, 0.29}
};

PolynomialPotential get_potential(int order) {
    using namespace casadi;
    assert(2 <= order <= 8);

    SXVector phi;
    SX delta = SX::sym("delta");

    std::vector<double> order_coeffs = coeffs[order - 2];

    for (int i = 0; i < order; ++i) {
        std::ostringstream name;
        name << "phi_" << order;
        phi.push_back(SX::sym(name.str()));
    }

    SX V_term_1 = 0;
    SX V_term_2 = 0;

    for (int i = 0; i < order; ++i) {
        V_term_1 += order_coeffs[i]*sqr(phi[i] - 1);
    }
    V_term_1 -= delta;

    for (int i = 0; i < order; ++i) {
        V_term_2 += sqr(phi[i]);
    }

    SX V = V_term_1*V_term_2;
    Function fV = Function("fV", {SX::vertcat(phi), delta}, {V}, {"phi", "delta"}, {"V"});
    
    PolynomialPotential pp;
    pp.fV = fV;
    pp.phi = SX::vertcat(phi);
    pp.delta = delta;
    return pp;
}

std::vector<double> find_true_vac(PolynomialPotential pp, double delta) {
    using namespace casadi;
    
    int n_phi = pp.phi.size1();
    DM ubarg = DM::vertcat(DMVector(n_phi, 3.));
    DM lbarg = DM::vertcat(DMVector(n_phi, 0.5));
    DM arg0 = DM::vertcat(DMVector(n_phi, 1.));

    return BubbleProfiler2::find_minimum(
        pp.fV, pp.phi, lbarg, ubarg, arg0, {pp.delta}, {delta}
    );
}

#endif