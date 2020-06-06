#include <BubbleProfiler2/CasadiBounceSolver.hpp>
#include <BubbleProfiler2/BouncePath.hpp>
#include <BubbleProfiler2/CasadiFindMinimum.hpp>

#define sqr(x) (x)*(x)

int main() {
    using namespace casadi;
    using namespace BubbleProfiler2;

    SX phi_1 = SX::sym("phi_1");
    SX phi_2 = SX::sym("phi_2");
    SX phi = SX::vertcat({phi_1, phi_2});

    SX delta = SX::sym("delta");
    DMVector param_vals = {0.4};

    SX V = (sqr(phi_1) + sqr(phi_2))*(1.8*sqr(phi_1 - 1) + 0.2*sqr(phi_2 - 1) - delta);
    
    Function fV = Function("fV", {phi, delta}, {V}, {"phi", "delta"}, {"V"});
    
    DM lbarg = DM::vertcat({0.5, 0.5});
    DM ubarg = DM::vertcat({2., 2.});
    DM arg0 = DM::vertcat({1., 1.});
    std::vector<double> true_vac = find_minimum(fV, phi, lbarg, ubarg, arg0, {delta}, param_vals);
    std::cout << true_vac << std::endl;
}
