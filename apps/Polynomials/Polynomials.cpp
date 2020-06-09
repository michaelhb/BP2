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
    {0.2434, 0.5233, 0.34234, 0.4747, 0.234808, 0.57023, 0.138912, 0.51723}
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

int main() {
    using namespace BubbleProfiler2;
    using namespace casadi;

    // Degree of the polynomial potential
    int order = 5;

    // Thinness parameter (smaller values -> thinner walled)
    double delta = 0.8;

    PolynomialPotential pp = get_potential(order);
    std::vector<double> true_vac = find_true_vac(pp, delta);
    std::vector<double> false_vac(order, 0.);

    CasadiBounceSolver solver = CasadiBounceSolver(pp.fV, order, {pp.delta}, 3, 100);

    std::map<std::string, double> v_pars;
    v_pars["delta"] = delta;

    BouncePath path = solver.solve(true_vac, false_vac, v_pars);
    
    std::cout << "Bounce action: " << path.get_action() << std::endl;
    path.plot_profiles(20., "title");
}