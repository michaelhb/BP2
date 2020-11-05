#ifndef BUBBLEPROFILER2_BOUNCEPATH_HPP_INCLUDED
#define BUBBLEPROFILER2_BOUNCEPATH_HPP_INCLUDED

#include <cassert>
#include <vector>
#include <string>
#include <map>

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
        metadata = std::map<std::string,double>();
    };

    // Just expose the data for now, add more appropriate 
    // accessors later as required.

    const std::vector<double>& get_radii() const {return radii;}

    const std::vector<std::vector<double>>& get_profiles() const {return profiles;}

    double get_action() const {return action;}

    void plot_profiles(double r_max = -1., std::string title = "") const;

    const double get_metadata(std::string label) const {
        return metadata.at(label);
    }

    void set_metadata(std::string label, double value) {
        metadata[label] = value;
    }

private:
    double action;
    std::vector<double> radii;
    std::vector<std::vector<double>> profiles;
    std::map<std::string, double> metadata;
};

};

#endif