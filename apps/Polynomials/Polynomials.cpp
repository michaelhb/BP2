
#include "Polynomials.hpp"

int main() {
    using namespace BubbleProfiler2;
    using namespace casadi;

    // Degree of the polynomial potential
    int order = 2;

    // Thinness parameter (smaller values -> thinner walled)
    double delta = 0.1;

    PolynomialPotential pp = get_potential(order);
    std::vector<double> true_vac = find_true_vac(pp, delta);
    std::vector<double> false_vac(order, 0.);

    CasadiBounceSolver solver = CasadiBounceSolver(pp.fV, order, {pp.delta}, 3, 100, false);

    std::map<std::string, double> v_pars;
    v_pars["delta"] = delta;

    BouncePath path = solver.solve(true_vac, false_vac, v_pars);
    
    std::cout << "Bounce action: " << path.get_action() << std::endl;
    path.plot_profiles(20., "title");
}