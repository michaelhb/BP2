#ifndef BUBBLEPROFILER2_FINDMINIMUM_HPP_INCLUDED
#define BUBBLEPROFILER2_FINDMINIMUM_HPP_INCLUDED

#include <casadi/casadi.hpp>

namespace BubbleProfiler2 {

std::vector<double> find_minimum(casadi::Function potential, 
    casadi::SX arg, casadi::DM lbarg, casadi::DM ubarg, casadi::DM arg0,
    casadi::SXVector params, casadi::DMVector param_vals) {
    using namespace casadi;
    
    // Create the NLP
    SXDict argV;
    argV["phi"] = arg;
    for (int i = 0; i < params.size(); ++i) {
        argV[params[i].name()] = params[i];
    }

    SX f = potential(argV).at("V");        
    SXDict nlp_arg = {{"f", f}, {"x", arg}, {"p", vertcat(params)}};
    Dict nlp_opt = Dict();
    nlp_opt["ipopt.print_level"] = 0;
    nlp_opt["print_time"] = 0;
    nlp_opt["ipopt.sb"] = "yes";
    Function solver = nlpsol("nlpsol", "ipopt", nlp_arg, nlp_opt);
    
    // Solve it
    DMDict sol_arg = {{"x0", arg0}, {"lbx", lbarg}, {"ubx", ubarg}, {"p", param_vals}};
    DMDict res = solver(sol_arg);
    return res.at("x").get_elements();
}

}

#endif