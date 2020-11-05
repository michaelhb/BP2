#include "Polynomials.hpp"
#include <chrono>

int main(int argc, char* argv[]) {
    using namespace BubbleProfiler2;
    using namespace std::chrono;

    int n_grid = 50;

    std::vector<double> deltas = {0.065, 0.11, 0.13, 0.15, 0.2, 0.22, 0.29, 0.27,
         0.3, 0.32, 0.39, 0.39, 0.42, 0.45, 0.52, 0.47, 0.56, 0.55};

    std::map<std::string, double> v_pars;

    for (int i = 0; i < deltas.size(); ++i) {
        int n_fields = i + 3;
        std::vector<double> false_vac(n_fields, 0.);
        double delta = deltas[i];
        // double delta = 0.15;
        v_pars["delta"] = delta;

        PolynomialPotential pp = get_potential(n_fields);
        auto t_setup_start = high_resolution_clock::now();
        CasadiBounceSolver solver = CasadiBounceSolver(pp.fV, n_fields, {pp.delta}, 3, n_grid, true);
        auto t_setup_end = high_resolution_clock::now();
        auto setup_duration = duration_cast<microseconds>(t_setup_end - t_setup_start).count() * 1e-6;

        std::vector<double> true_vac = find_true_vac(pp, delta);
        BouncePath path = solver.solve(true_vac, false_vac, v_pars);

        double solve_duration = path.get_metadata("solve_duration");
        
        std::cout << n_fields << "," << delta << "," << path.get_action() << "," << solve_duration << "," << setup_duration << std::endl;
    }

}