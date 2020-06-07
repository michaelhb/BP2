#include <BubbleProfiler2/BouncePath.hpp>
#include "gnuplot-iostream.h"

namespace BubbleProfiler2 {

void BouncePath::plot_profiles(double r_max, std::string title) const {
        int n_phi = profiles[0].size();
        int n_rho = profiles.size();

        std::vector<std::vector<double>> profile_data;

        for (int i = 0; i < n_rho; ++i) {
            if (r_max > 0 && radii[i] > r_max) {
                break;
            }
            std::vector<double> row = {radii[i]};
            for (int j = 0; j < n_phi; ++j) {
                row.push_back(profiles[i][j]);
            }
            profile_data.push_back(row);
        }        

        Gnuplot gp;
        gp << "set title '" << title << "'\n";
        gp << "plot '-' using 1:2 with lines title 'phi_1'";
        for (int i = 2; i <= n_phi; ++i) {
            gp << ", '' using 1:" << i + 1 << " with lines title 'phi_" << i << "'";
        }
        gp << "\n";

        for (int i = 0; i < n_phi; ++i) {
            gp.send1d(profile_data);
        }
}

}