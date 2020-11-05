#include <BubbleProfiler2/CasadiBounceSolver.hpp>

namespace BubbleProfiler2 {

CasadiBounceSolver::CasadiBounceSolver(casadi::Function potential_, int n_phi_,
        casadi::SXVector v_params_, int n_dims_, int N_, bool quiet_) :
        n_phi(n_phi_), n_dims(n_dims_), N(N_), v_params(v_params_), potential(potential_), quiet(quiet_) {
        using namespace casadi;

        // TODO most of what's below probably shouldn't live in the 
        // constructor...

        // Get the internal coordinate grid for each finite element
        tau_root = collocation_points(d, "legendre");
        tau_root.insert(tau_root.begin(), 0.);

        // Control intervals and element sizes (evenly spaced for now)
        for (int i = 0; i < N; ++i) {
            double t = -1.0 + 2.0*i / N;
            double h = 2.0/N;
            t_k.push_back(t);
            h_k.push_back(h);
        }

        // Surface area of N sphere (TODO look up general formula?)
        if (n_dims == 3) {
            S_n = 4*pi;
        }
        else if (n_dims == 4) {
            S_n = 0.5*pi*pi;
        }
        else {
            throw std::invalid_argument("Only d = 3 and d = 4 are currently supported.");
        }

        // Initialise polynomial basis and collocation / integration coefficients
        // TODO - this should be done offline, since (for now) we stick with cubics
        D = std::vector<double>(d + 1, 0);
        B = std::vector<double>(d + 1, 0);
        P = std::vector<Polynomial>(d + 1);
        C = std::vector<std::vector<double> >(d+1,std::vector<double>(d+1, 0));

        // Construct polynomial basis & extract relevant coefficients
        for (int j = 0; j < d + 1; ++j) {

            Polynomial p = 1;
            for(int r = 0; r < d + 1; ++r){
                if(r != j){
                    p *= Polynomial(-tau_root[r],1)/(tau_root[j]-tau_root[r]);
                }
            }
            P[j] = p;

            // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
            D[j] = p(1.0);

            // Evaluate the time derivative of the polynomial at all collocation points 
            Polynomial dp = p.derivative();
            for(int r=0; r<d+1; ++r){
                C[j][r] = dp(tau_root[r]);
            }
            
            // Evaluate the integral of the polynomial to get the coefficients for Gaussian quadrature
            Polynomial ip = p.anti_derivative();
            B[j] = ip(1.0);
        }

        // Construct the NLP
        nlp = get_nlp(potential);
    }

Static_bounds CasadiBounceSolver::get_static_bounds(std::vector<double> false_vac) const {
    using namespace casadi;
    /**** Bounds on control variables ****/
    // Need to find a less opaque way of ensuring that 
    // concatenated NLP inputs / bounds / start values 
    // are consistently ordered!

    // Limits for unbounded variables
    std::vector<double> ubinf(n_phi, inf);
    std::vector<double> lbinf(n_phi, -inf);

    // Zero vector for constraint bounds
    std::vector<double> zeroes(n_phi, 0);
    
    // Zero vector for collocation bounds
    std::vector<double> zeroes_col(d*n_phi, 0);
    std::vector<double> lbU, ubU;

    // Derivative at origin fixed to zero
    append_d(lbU, zeroes);
    append_d(ubU, zeroes);

    for (int k = 1; k <= N; ++k) {
        append_d(lbU, lbinf);
        append_d(ubU, ubinf);
    }

    /**** Bounds on state variables ****/
    std::vector<double> lbPhi, ubPhi;

    // Free endpoint states
    for (int k = 0; k < N; ++k) {
        append_d(lbPhi, lbinf);
        append_d(ubPhi, ubinf);
    }

    // Final state, fixed to the false vacuum
    append_d(lbPhi, false_vac);
    append_d(ubPhi, false_vac);

    // Free intermediate states
    for (int k = 0; k < N; ++k) {
        for (int j = 1; j <= d; ++j) {
            append_d(lbPhi, lbinf);
            append_d(ubPhi, ubinf);
        }
    }

    /**** Bounds on constraints ****/
    std::vector<double> lbg = {}; // Lower bounds for constraints
    std::vector<double> ubg = {}; // Upper bounds for constraints

    // Continuity equations
    for (int k = 0; k < N; ++k) {
        append_d(lbg, zeroes);
        append_d(ubg, zeroes);
    }

    // Collocation equations and objective function
    for (int k = 0; k < N; ++k) {
        append_d(lbg, zeroes_col);
        append_d(ubg, zeroes_col);
    }

    // Add potential constraint
    lbg.push_back(0);
    ubg.push_back(0);

    // Return results
    Static_bounds bounds;
    bounds.lbU = lbU;
    bounds.ubU = ubU;
    bounds.lbPhi = lbPhi;
    bounds.ubPhi = ubPhi;  
    bounds.lbg = lbg;
    bounds.ubg = ubg;

    return bounds;
}

NLP CasadiBounceSolver::get_nlp(const casadi::Function& potential) const {
    using namespace casadi;
    using namespace std::chrono;

    auto t_setup_start = high_resolution_clock::now();

    /**** Initialise parameter variables ****/
    SXVector h_par, gamma_par, gammadot_par; 
    SXVector grid_pars; // All grid parameter variables

    for (int i = 0; i < N; ++i) {
        SX h_par_ = SX::sym(varname("h", {i}));
        h_par.push_back(h_par_);
        grid_pars.push_back(h_par_);
    }
    for (int i = 0; i < N; ++i) {
        SX gamma_par_ = SX::sym(varname("gamma", {i}), d + 1);
        gamma_par.push_back(gamma_par_);
        grid_pars.push_back(gamma_par_);
    }
    for (int i = 0; i < N; ++i) {
        SX gammadot_par_ = SX::sym(varname("gammadot", {i}), d + 1);
        gammadot_par.push_back(gammadot_par_);
        grid_pars.push_back(gammadot_par_);
    }

    SX V0_par = SX::sym("V0");
    
    /**** Symbolic ansatz parameters ****/
    SX r0 = SX::sym("r0");
    SX sigma = SX::sym("sigma");
    SX true_vac = SX::sym("true_vacuum", n_phi);
    SX false_vac = SX::sym("false_vacuum", n_phi);
    SX zero_field = SX::sym("zero", n_phi);

    /**** Initialise control variables ****/        
    SXVector U, U_ansatz;

    // Derivative at origin fixed to zero
    SX U_0_0 = SX::sym("U_0_0", n_phi);
    U.push_back(U_0_0);
    U_ansatz.push_back(ansatz_dot(0, 0, true_vac, false_vac, r0, sigma, gamma_par));

    for (int k = 1; k < N; ++k) {
        SX Uk = SX::sym(varname("U", {k}), n_phi);
        U.push_back(Uk);
        U_ansatz.push_back(ansatz_dot(k, 0, true_vac, false_vac, r0, sigma, gamma_par));
    }

    // Endpoints (avoid Tr(1) singularity)
    SX UN = SX::sym(varname("U", {N}), n_phi);
    U.push_back(UN);
    U_ansatz.push_back(zero_field); // Does this make any sense?!

    /**** Initialise state variables ****/
    SXVector Phi, Phi_ansatz;
    std::vector<SXVector> element_states; // States within an element
    SXVector element_plot; // Concatenated states within an element
    SXVector endpoints; // Endpoint states

    // Free endpoint states
    for (int k = 0; k < N; ++k) {
        SX phi_k_0 = SX::sym(varname("phi", {k, 0}), n_phi);
        endpoints.push_back(phi_k_0);
        Phi.push_back(phi_k_0);
        Phi_ansatz.push_back(ansatz(k, 0, true_vac, false_vac, r0, sigma, gamma_par));
    }

    // Final state, fixed to the false vacuum
    SX phi_N_0 = SX::sym("phi_N_0", n_phi);
    endpoints.push_back(phi_N_0);
    Phi.push_back(phi_N_0);
    Phi_ansatz.push_back(false_vac);

    // Build finite elements (including left endpoints)
    for (int k = 0; k < N; ++k) {
        std::vector<SX> e_states;
        e_states.push_back(endpoints[k]);
        for (int j = 1; j <= d; ++j) {
            SX phi_k_j = SX::sym(varname("phi", {k, j}), n_phi);
            e_states.push_back(phi_k_j);
            Phi.push_back(phi_k_j);
            Phi_ansatz.push_back(ansatz(k, j, true_vac, false_vac, r0, sigma, gamma_par));
        }
        element_plot.push_back(SX::horzcat(e_states));
        element_states.push_back(e_states);
    }

    /**** Useful functions of the state and control variables in an element ****/

    // State variables in a given element
    SXVector element;
    for (int j = 0; j <= d; ++j) {
        element.push_back(SX::sym(varname("phi_k", {j}), n_phi));
    }        

    // Width of control interval
    SX h_elem = SX::sym("h_elem");

    // Monotonic transformation & derivative
    SX gamma = SX::sym("gamma", d + 1);
    SX gammadot = SX::sym("gammadot", d + 1);
    
    // Estimate for the state at end of control interval
    SX phi_end = 0;

    for (int i = 0; i <= d; ++i) {
        phi_end += D[i]*element[i];
    }
    
    Function Phi_end = Function("Phi_end", element, {phi_end});

    // Interpolated controls in an element
    SX control_start = SX::sym("u_k", n_phi);
    SX control_end = SX::sym("u_k+1", n_phi);
    SXVector control_int;

    for (int j = 1; j <= d; ++j) {
        control_int.push_back(
            (1 - tau_root[j])*control_start + tau_root[j]*control_end
        );
    }
    
    // Derivative constraints in an element
    SXVector phidot_cons;
    
    for (int j = 1; j <= d; ++j) {
        SX phidot_approx = 0;
        for (int r = 0; r <= d; ++r) {
            phidot_approx += C[r][j]*element[r];
        }
        phidot_cons.push_back(h_elem*gammadot(j)*control_int[j - 1] - phidot_approx);
    }

    SXVector phidot_inputs;
    phidot_inputs.insert(phidot_inputs.end(), element.begin(), element.end());
    phidot_inputs.push_back(control_start);
    phidot_inputs.push_back(control_end);
    phidot_inputs.push_back(h_elem);
    phidot_inputs.push_back(gammadot);
    
    Function Phidot_cons = Function("Phidot_cons", phidot_inputs, phidot_cons);

    // Value of kinetic objective on an element
    SX T_k = 0;

    for (int j = 1; j <= d; ++j) {
        T_k = T_k + 0.5*S_n*h_elem*B[j]*pow(gamma(j), n_dims - 1)
            *gammadot(j)*dot(control_int[j - 1], control_int[j - 1]);
    }

    // Value of potential constraint functional on an element
    SX V_k = 0;

    for (int j = 1; j <= d; ++j) {
        SXDict argV;
        argV["phi"] = element[j];
        add_param_args(argV);
        V_k = V_k + S_n*h_elem*B[j]*pow(gamma(j), n_dims - 1)*gammadot(j)*potential(argV).at("V");
    }

    // TODO clean up function inputs, this is messy / verbose

    // We define per-element functions for the quadratures, then use them to build
    // integrals over the whole domain.
    SXVector quadrature_inputs;
    quadrature_inputs.insert(quadrature_inputs.end(), element.begin(), element.end());
    quadrature_inputs.push_back(control_start);
    quadrature_inputs.push_back(control_end);
    quadrature_inputs.push_back(h_elem);
    quadrature_inputs.push_back(gamma);
    quadrature_inputs.push_back(gammadot);
    quadrature_inputs.insert(quadrature_inputs.end(), v_params.begin(), v_params.end());
    
    Function T_obj_k = Function("T_obj_k", quadrature_inputs, {T_k});
    Function V_cons_k = Function("V_cons_k", quadrature_inputs, {V_k});

    // Build quadratures
    SX T = 0; // Objective function
    SX V = 0; // Integral constraint

    for (int k = 0; k < N; ++k) {
        SXVector quadrature_inputs_ = SXVector(element_states[k]);
        quadrature_inputs_.push_back(U[k]);
        quadrature_inputs_.push_back(U[k + 1]);
        quadrature_inputs_.push_back(h_par[k]);
        quadrature_inputs_.push_back(gamma_par[k]);
        quadrature_inputs_.push_back(gammadot_par[k]);
        quadrature_inputs_.insert(quadrature_inputs_.end(), v_params.begin(), v_params.end());

        T += T_obj_k(quadrature_inputs_)[0];
        V += V_cons_k(quadrature_inputs_)[0];
    }
    
    SXVector V_a_arg = {vertcat(Phi), vertcat(grid_pars)};
    StringVector V_a_argnames = {"Phi", "par"};
    SXVector T_a_arg = {vertcat(U), vertcat(grid_pars)};
    StringVector T_a_argnames = {"U", "par"};

    for (int i = 0; i < v_params.size(); ++i) {
        V_a_arg.push_back(v_params[i]);
        V_a_argnames.push_back(v_params[i].name());
        T_a_arg.push_back(v_params[i]);
        T_a_argnames.push_back(v_params[i].name());
    }

    Function V_a = Function("fV", V_a_arg, {V}, V_a_argnames, {"V"});
    Function T_a = Function("fT", T_a_arg, {T}, T_a_argnames, {"T"});

    // We want these as functions of {r0,sigma,v_params,grid_pars}
    SXDict V_ansatz_args, T_ansatz_args;
    V_ansatz_args["Phi"] = vertcat(Phi_ansatz);
    V_ansatz_args["par"] = vertcat(grid_pars);
    T_ansatz_args["U"] = vertcat(U_ansatz);
    T_ansatz_args["par"] = vertcat(grid_pars);
    add_param_args(V_ansatz_args);
    add_param_args(T_ansatz_args);

    /**** Build constraints ****/ 
    SXVector g = {}; // All constraints

    // Continuity equations
    for (int k = 0; k < N; ++k) {
        g.push_back(Phi_end(element_states[k])[0] - endpoints[k + 1]);
    }

    // Collocation equations
    for (int k = 0; k < N; ++k) {
        SXVector phidot_inputs_ = SXVector(element_states[k]);
        phidot_inputs_.push_back(U[k]);
        phidot_inputs_.push_back(U[k + 1]);
        phidot_inputs_.push_back(h_par[k]);
        phidot_inputs_.push_back(gammadot_par[k]);
        g.push_back(SX::vertcat(Phidot_cons(phidot_inputs_)));
    }

    // Potential constraint
    g.push_back(V0_par - V);

    /**** Concatenate variables ****/
    SXVector w = {}; // All decision variables
    w.insert(w.end(), Phi.begin(), Phi.end());
    w.insert(w.end(), U.begin(), U.end());

    // Collect states and constraints into single vectors
    SX W = SX::vertcat(w);
    SX G = SX::vertcat(g);

    // Versions of T and V suitable for evaluating on results 
    SX elements_plot = SX::horzcat(element_plot);
    Function Phi_ret = Function("elements", {W}, {elements_plot});
    Function T_ret = Function("T_ret", {W, vertcat(grid_pars), vertcat(v_params)}, {T}, {"W", "Par", "v_params"}, {"T"});
    Function V_ret = Function("V_ret", {W, vertcat(grid_pars), vertcat(v_params)}, {V}, {"W", "Par", "v_params"}, {"V"});

    /**** Create the main NLP ****/
    
    // Grid pars + V0 + v_pars
    SXVector pars = grid_pars;
    pars.push_back(V0_par);
    for (auto & v_param: v_params) {
        pars.push_back(v_param);
    }

    SXDict nlp_arg = {{"f", T}, {"x", W}, {"g", G}, {"p", vertcat(pars)}};
    Dict nlp_opt = Dict();
    nlp_opt["expand"] = false;
    nlp_opt["ipopt.linear_solver"] = "ma27";
    
    if (quiet) {
        nlp_opt["ipopt.print_level"] = 0;
        nlp_opt["print_time"] = 0;
        nlp_opt["ipopt.sb"] = "yes";
    }

    Function solver = nlpsol("full_nlp", "ipopt", nlp_arg, nlp_opt);

    /**** Create the ansatz NLP ****/
    SXVector ansatz_nlp_pars(pars);
    ansatz_nlp_pars.push_back(true_vac);
    ansatz_nlp_pars.push_back(false_vac);
    ansatz_nlp_pars.push_back(zero_field); // Hacky :(

    Dict ansatz_nlp_opt(nlp_opt);

    SXDict ansatz_nlp_arg = {
        {"f", T_a(T_ansatz_args).at("T")}, 
        {"x", SX::vertcat({r0, sigma})},
        {"g", V0_par - V_a(V_ansatz_args).at("V")},
        {"p", vertcat(ansatz_nlp_pars)}};

    Function ansatz_solver = nlpsol("ansatz_nlp", "ipopt", ansatz_nlp_arg, ansatz_nlp_opt);

    auto t_setup_end = high_resolution_clock::now();
    auto setup_duration = duration_cast<microseconds>(t_setup_end - t_setup_start).count() * 1e-6;
    if (!quiet) std::cout << "CasadiBounceSolver - setup took " << setup_duration << " sec" << std::endl;
    
    NLP nlp;
    nlp.ansatz_nlp = ansatz_solver;
    nlp.nlp = solver;
    nlp.T_a = T_a;
    nlp.T_ret = T_ret;
    nlp.V_a = V_a;
    nlp.V_ret = V_ret;
    nlp.Phi_ret = Phi_ret;

    return nlp;
}

Ansatz CasadiBounceSolver::get_initial_ansatz(
    casadi::Function fV, 
    casadi::Function fT,
    CompactGrid& grid,
    std::map<std::string, double> v_pars, 
    casadi::DM true_vac, casadi::DM false_vac) const {

    Ansatz a;
    std::vector<double> Phi0, U0;
    Phi0.reserve((N*(n_dims + 1) + 1));
    U0.reserve(N + 1);
    casadi::DMDict argV(v_pars.begin(), v_pars.end());
    argV["par"] = grid.concatenate(); 

    // CURRENT BEST SETTINGS 
    double r0 = 0.5;
    double delta_r0 = 0.5;
    double r0_max = 25.;
    double targetV = -1;
    double tol = 1e-3;
    double sig_upper0 = 3.;
    double sig_upper = sig_upper0;
    double sig_lower = 0.;
    double sigma_min = 0.05;
    double sigma = 0.5*(sig_upper + sig_lower);
    double sigma_cache = sigma;
    double V_mid;

    double scale = grid.grid_scale();

    do {
        Phi0.clear();

        // Endpoints
        for (int k = 0; k <= N; ++k) {
            append_d(Phi0, ansatz(
                Tr(t_kj(k, 0), scale), true_vac, false_vac, r0, sigma).get_elements());
        }

        // Collocation points
        for (int k = 0; k < N; ++k) {
            for (int j = 1; j <= d; ++j) {
                append_d(Phi0, ansatz(
                    Tr(t_kj(k, j), scale), true_vac, false_vac, r0, sigma ).get_elements());
            }
        }

        argV["Phi"] = Phi0;
        
        V_mid = fV(argV).at("V").get_elements()[0];

        if (!quiet) {
            std::cout << "r0 = " << r0 << ", sigma = " << sigma 
                        << ", V_mid = " << V_mid << std::endl; 
        }

        if (V_mid < targetV) {
            sig_lower  = sigma;
        }
        else if (V_mid > targetV) {
            sig_upper = sigma;
        }
        sigma_cache = sigma;
        sigma = 0.5*(sig_lower + sig_upper);
        
        if (sigma < sigma_min) {
            r0 += delta_r0;
            sig_upper = sig_upper0;
            sig_lower = 0;
            sigma = 0.5*(sig_lower + sig_upper);
            V_mid = 1; // Make sure we loop again

            if (r0 > r0_max) {
                throw std::runtime_error("Exceeded maximum radius!");
            }
        }
    }
    while (std::abs(V_mid - targetV) > tol);

    // Calculate the derivatives
    for (int k = 0; k < N; ++k) {
        append_d(U0, ansatz_dot(
            Tr(t_kj(k,0), scale), true_vac, false_vac, r0, sigma_cache).get_elements());
    } 

    // Avoid Tr(1) singularity
    append_d(U0, std::vector<double>(false_vac.size1(), 0));

    if (!quiet) {
        std::cout << "Ansatz r0 = " << r0 << ", sigma = " << sigma_cache << std::endl;
    }

    a.V0 = V_mid;
    a.r0 = r0;
    a.sigma = sigma_cache;
    a.true_vac = true_vac;
    a.false_vac = false_vac;
    a.Phi0 = Phi0;
    a.U0 = U0;
    return a;
}

Ansatz CasadiBounceSolver::rescale_ansatz(Ansatz a, CompactGrid& grid) const {
    std::vector<double> Phi0, U0;
    double scale = grid.grid_scale();

    // State
    for (int k = 0; k <= N; ++k) {
        append_d(Phi0, ansatz(
            Tr(t_kj(k, 0), scale), a.true_vac, a.false_vac, a.r0, a.sigma).get_elements());
    }

    // State collocation points
    for (int k = 0; k < N; ++k) {
        for (int j = 1; j <= d; ++j) {
            append_d(Phi0, ansatz(
                Tr(t_kj(k, j), scale), a.true_vac, a.false_vac, a.r0, a.sigma).get_elements());
        }
    }

    // Calculate the derivatives
    for (int k = 0; k < N; ++k) {
        append_d(U0, ansatz_dot(
            Tr(t_kj(k,0), scale), a.true_vac, a.false_vac, a.r0, a.sigma).get_elements());
    } 

    // Avoid Tr(1) singularity
    append_d(U0, std::vector<double>(a.false_vac.size1(), 0));    

    Ansatz a_scaled;
    a_scaled.V0 = a.V0;
    a_scaled.r0 = a.r0;
    a_scaled.sigma = a.sigma;
    a_scaled.true_vac = a.true_vac;
    a_scaled.false_vac = a.false_vac;
    a_scaled.Phi0 = Phi0;
    a_scaled.U0 = U0;

    return a_scaled;
}

Ansatz CasadiBounceSolver::get_ansatz(double r0, double sigma, double V0, std::vector<double> true_vac, 
        std::vector<double> false_vac, CompactGrid grid) const {
    // TODO: Make Ansatz a class already, this is terrible
    Ansatz a;
    a.r0 = r0;
    a.sigma = sigma;
    a.true_vac = true_vac;
    a.false_vac = false_vac;
    a.V0 = V0;
    
    return rescale_ansatz(a, grid);
}

BouncePath CasadiBounceSolver::solve(const std::vector<double>& true_vacuum, const std::vector<double>& false_vacuum,
    std::map<std::string, double> v_pars) const {

    using namespace casadi;
    using namespace std::chrono;

    if (!quiet) std::cout << "TRUE_VACUUM: " << true_vacuum << std::endl;
    if (!quiet) std::cout << "FALSE_VACUUM: " << false_vacuum << std::endl;

    // Get concatenated grid parameters (h_k, gamma, gamma_dot)
    CompactGrid grid = get_grid(default_grid_scale);

    /**** Find the initial ansatz ****/
    auto t_ansatz_start = high_resolution_clock::now();
    Ansatz initial_ansatz = get_initial_ansatz(nlp.V_a, nlp.T_a, grid, v_pars, true_vacuum, false_vacuum);

    if (!quiet) std::cout << "Initial ansatz r0 = " << initial_ansatz.r0 
        << ", sigma = " << initial_ansatz.sigma << std::endl;

    // Evaluate T and V on the ansatz
    DMDict argT(v_pars.begin(), v_pars.end());
    argT["U"] = initial_ansatz.U0;
    argT["par"] = grid.concatenate();
    double T0 = nlp.T_a(argT).at("T").get_elements()[0];

    DMDict argV(v_pars.begin(), v_pars.end());
    argV["Phi"] = initial_ansatz.Phi0;
    argV["par"] = grid.concatenate(); 
    double V0 = nlp.V_a(argV).at("V").get_elements()[0];

    // Create an adjusted grid to increase the number of points 
    // near the expected bubble wall location
    double new_scale = find_grid_scale(initial_ansatz.r0);

    grid = get_grid(new_scale);
    if (!quiet) std::cout << "Adjusted grid scale to " << new_scale << std::endl;

    // Rescale the ansatz to the new grid
    Ansatz scaled_ansatz = rescale_ansatz(initial_ansatz, grid);

    /**** Refine the ansatz ****/

    // Vector of (canonically ordered) model parameter values
    std::vector<double> v_param_vals;
    for (auto & v_param_sx : v_params) {
        v_param_vals.push_back(v_pars[v_param_sx.name()]);
    }

    // Concat the grid params, V0 param value, and model params
    std::vector<double> pars(grid.concatenate());
    pars.push_back(V0);
    pars.insert(pars.end(), v_param_vals.begin(), v_param_vals.end());

    // pars will be used for the full NLP, but we need some extra 
    // params for the ansatz
    std::vector<double> ansatz_pars(pars);
    ansatz_pars.insert(ansatz_pars.end(), true_vacuum.begin(), true_vacuum.end());
    ansatz_pars.insert(ansatz_pars.end(), false_vacuum.begin(), false_vacuum.end());
    for (int i = 0; i < n_phi; ++i) ansatz_pars.push_back(0); // Hacky :(

    std::vector<double> x0_a = {scaled_ansatz.r0, scaled_ansatz.sigma};
    // std::vector<double> lbx_a = {1e-5, 1e-5};
    std::vector<double> lbx_a = {1e-3, 0.02};
    std::vector<double> ubx_a = {100., 100.};
    std::vector<double> lbg_a = {0.};
    std::vector<double> ubg_a = {0.};

    DMDict ansatz_arg = {{"x0", x0_a}, {"lbx", lbx_a}, {"ubx", ubx_a}, 
        {"lbg", lbg_a}, {"ubg", ubg_a}, {"p", ansatz_pars}};

    DMDict res_a = nlp.ansatz_nlp(ansatz_arg);
    double r0_a = res_a.at("x").get_elements()[0];
    double sigma_a = res_a.at("x").get_elements()[1];

    if (!quiet) std::cout << "Final ansatz r0 = " << r0_a 
        << ", sigma = " << sigma_a << std::endl;

    Ansatz final_ansatz = get_ansatz(r0_a, sigma_a, V0, true_vacuum, false_vacuum, grid);
    // Ansatz final_ansatz = scaled_ansatz;

    auto t_ansatz_end = high_resolution_clock::now();
    auto ansatz_duration = duration_cast<microseconds>(t_ansatz_end - t_ansatz_start).count() * 1e-6;
    if (!quiet) std::cout << "Finding ansatz took " << ansatz_duration << "s" << std::endl;

    /**** Concatenate NLP initial state ****/
    std::vector<double> w0 = {}; // Initial values for decision variables
    w0.insert(w0.end(), final_ansatz.Phi0.begin(), final_ansatz.Phi0.end());
    w0.insert(w0.end(), final_ansatz.U0.begin(), final_ansatz.U0.end());

    /**** Concatenate decision variable bounds****/
    Static_bounds bounds = get_static_bounds(false_vacuum);
    std::vector<double> lbw = {}; 
    std::vector<double> ubw = {};  

    lbw.insert(lbw.end(), bounds.lbPhi.begin(), bounds.lbPhi.end());
    ubw.insert(ubw.end(), bounds.ubPhi.begin(), bounds.ubPhi.end());
    lbw.insert(lbw.end(), bounds.lbU.begin(), bounds.lbU.end());
    ubw.insert(ubw.end(), bounds.ubU.begin(), bounds.ubU.end());

    /**** Initialise and solve the NLP ****/

    // Run the optimiser. This is the other bottleneck, so we time it too.
    DMDict arg = {{"x0", w0}, {"lbx", lbw}, {"ubx", ubw}, {"lbg", bounds.lbg}, {"ubg", bounds.ubg}, {"p", pars}};
    auto t_solve_start = high_resolution_clock::now();
    DMDict res = nlp.nlp(arg);
    auto t_solve_end = high_resolution_clock::now();
    auto solve_duration = duration_cast<microseconds>(t_solve_end - t_solve_start).count() * 1e-6;
    if (!quiet) std::cout << "CasadiMaupertuisSolver - optimisation took " << solve_duration << " sec" << std::endl;

    // Evaluate the objective & constraint on the result
    DMDict ret_arg;
    ret_arg["W"] = res["x"];
    ret_arg["Par"] = grid.concatenate();
    ret_arg["v_params"] = v_param_vals;

    double Tret = nlp.T_ret(ret_arg).at("T").get_elements()[0];
    double Vret = nlp.V_ret(ret_arg).at("V").get_elements()[0];

    if (!quiet) std::cout << "V(result) = " << Vret << std::endl;
    if (!quiet) std::cout << "T(result) = " << Tret << std::endl;
    
    // Calculate the action
    double action = std::pow(((2.0 - d)/d)*(Tret/Vret), 0.5*d)*((2.0*Vret)/(2.0 - d));

    // Return the result
    std::vector<double> radii(N*(d + 1));
    int c = 0;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j <= d; ++j) {
            radii[c] = Tr(t_kj(k, j), grid.grid_scale());
            c++;
        }
    }

    std::vector<double> phi_ret = nlp.Phi_ret(res["x"]).at(0).get_elements();
    std::vector<std::vector<double>> profiles(radii.size());

    c = 0;
    for (int row = 0; row < radii.size(); ++row) {
        profiles[row] = std::vector<double>(n_phi);
        for (int col = 0; col < n_phi; ++col) {
            profiles[row][col] = phi_ret[c];
            c++;
        }
    }

    BouncePath path = BouncePath(radii, profiles, action);
    path.set_metadata("solve_duration", solve_duration);
    return path;
}

};