#ifndef BUBBLEPROFILER2_BOUNCEPATH_HPP_INCLUDED
#define BUBBLEPROFILER2_BOUNCEPATH_HPP_INCLUDED

#include <cassert>
#include <vector>
#include <string>

namespace BubbleProfiler2 {

//! Class to represent results of the bounce calculation
class BouncePath {
public:
    BouncePath() = default;
    
    BouncePath(const std::vector<double> radii_,
               const std::vector<std::vector<double>> profiles_,
               double action_) {
        assert(radii_.size() == profiles_.size());
        radii = radii_;
        profiles = profiles_;
        action = action_;
    };

    // Just expose the data for now, add more appropriate 
    // accessors later as required.

    const std::vector<double>& get_radii() const {return radii;}

    const std::vector<std::vector<double>>& get_profiles() const {return profiles;}

    double get_action() const {return action;}

    void plot_profiles(double r_max = -1., std::string title = "") const;

private:
    double action;
    std::vector<double> radii;
    std::vector<std::vector<double>> profiles;
};

};

#endif