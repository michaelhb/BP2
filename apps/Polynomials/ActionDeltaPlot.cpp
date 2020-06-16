#include "Polynomials.hpp"
#include <chrono>

int main(int argc, char* argv[]) {
    using namespace BubbleProfiler2;
    using namespace std::chrono;
    
    // CLI params: number of fields, number of grid points, delta_min, delta_max, delta_step
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0] << "test_tag n_fields n_grid delta_min delta_max delta_step" << std::endl;
        return 1;
    }

    std::string test_tag = argv[1];
    int n_fields = atoi(argv[2]);
    int n_grid = atoi(argv[3]);
    double delta_min = atof(argv[4]);
    double delta_max = atof(argv[5]);
    double delta_step = atof(argv[6]);

    // Output format will be: test_tag,delta,S_E,time.
    // Not including header, since I'll want to concatenate these.

    std::vector<double> false_vac(n_fields, 0.);

    PolynomialPotential pp = get_potential(n_fields);
    CasadiBounceSolver solver = CasadiBounceSolver(pp.fV, n_fields, {pp.delta}, 3, n_grid, true);

    double delta = delta_min;
    std::map<std::string, double> v_pars;
    
    while (delta < delta_max) {
        v_pars["delta"] = delta;

        std::vector<double> true_vac = find_true_vac(pp, delta);

        auto t_solve_start = high_resolution_clock::now();
        BouncePath path = solver.solve(true_vac, false_vac, v_pars);
        auto t_solve_end = high_resolution_clock::now();
        auto solve_duration = duration_cast<microseconds>(t_solve_end - t_solve_start).count() * 1e-6;

        std::cout << '"' << test_tag << '"' << "," << delta << "," << path.get_action() << "," << solve_duration << std::endl;
        delta += delta_step;
    }
}