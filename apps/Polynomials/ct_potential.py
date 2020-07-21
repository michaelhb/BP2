import sys, os
import time
import numpy as np
from scipy import optimize
from collections import namedtuple
from cosmoTransitions import generic_potential, pathDeformation

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

coeffs = [
    [1.8, 0.2],
    [0.684373, 0.181928, 0.295089],
    [0.534808, 0.77023, 0.838912, 0.00517238],
    [0.4747, 0.234808, 0.57023, 0.138912, 0.517238],
    [0.34234, 0.4747, 0.234808, 0.57023, 0.138912, 0.517238],
    [0.5233, 0.34234, 0.4747, 0.234808, 0.57023, 0.138912, 0.517238],
    [0.2434, 0.5233, 0.34234, 0.4747, 0.234808, 0.57023, 0.138912, 0.51723],
    [0.21, 0.24, 0.52, 0.34, 0.47, 0.23, 0.57, 0.14, 0.52],
    [0.12, 0.21, 0.24, 0.52, 0.34, 0.47, 0.23, 0.57, 0.14, 0.52],
    [0.23, 0.21, 0.21, 0.24, 0.52, 0.34, 0.47, 0.23, 0.57, 0.14, 0.52],
    [0.12, 0.11, 0.12, 0.21, 0.24, 0.52, 0.34, 0.47, 0.23, 0.57, 0.14, 0.52],
    [0.54, 0.47, 0.53, 0.28, 0.35, 0.27, 0.42, 0.59, 0.33, 0.16, 0.38, 0.35, 0.17],
    [0.39, 0.23, 0.26, 0.40, 0.11, 0.42, 0.41, 0.27, 0.42, 0.54, 0.18, 0.59, 0.13, 0.29],
    [0.21, 0.22, 0.22, 0.23, 0.39, 0.55, 0.43, 0.12, 0.16, 0.58, 0.25, 0.50, 0.45, 0.35, 0.45],
    [0.42, 0.34, 0.43, 0.22, 0.59, 0.41, 0.58, 0.41, 0.26, 0.45, 0.16, 0.31, 0.39, 0.57, 0.43, 0.10],
    [0.24, 0.35, 0.39, 0.56, 0.37, 0.41, 0.52, 0.31, 0.52, 0.22, 0.58, 0.39, 0.39, 0.17, 0.46, 0.30, 0.37],
    [0.18, 0.17, 0.30, 0.22, 0.38, 0.48, 0.11, 0.49, 0.43, 0.47, 0.21, 0.29, 0.32, 0.36, 0.30, 0.56, 0.46, 0.42],
    [0.40, 0.14, 0.10, 0.43, 0.39, 0.27, 0.33, 0.59, 0.48, 0.36, 0.24, 0.28, 0.51, 0.59, 0.40, 0.39, 0.24, 0.35, 0.20],
    [0.42, 0.11, 0.47, 0.13, 0.16, 0.24, 0.58, 0.53, 0.38, 0.44, 0.18, 0.46, 0.47, 0.27, 0.53, 0.24, 0.33, 0.40, 0.32, 0.29]
]

class PolynomialPotential(generic_potential.generic_potential):

    def init(self, Ndim, delta):
        self.Ndim = Ndim
        self.delta = delta
        self.order_coeffs = coeffs[Ndim - 2]

    def V0(self, X):
        X = np.asanyarray(X)
        
        fields = []
        for i in range(self.Ndim):
            fields.append(X[...,i])

        v_term_1 = 0
        v_term_2 = 0
        
        for i in range(self.Ndim):
            v_term_1 += self.order_coeffs[i]*(fields[i] - 1)*(fields[i] - 1)
        
        v_term_1 -= self.delta

        for i in range(self.Ndim):
            v_term_2 += fields[i]*fields[i]

        return v_term_1*v_term_2

Point = namedtuple("Point", ["potential", "delta", "true_vac", "false_vac"])

def gen_points(Ndim, delta_min, delta_max, step):
    points = []
    delta = delta_min
    
    while delta < delta_max:
        potential = PolynomialPotential(Ndim, delta)
        v_ = lambda phi: potential.Vtot(phi,0)
        bnd = ((.5, 3),)*Ndim
        phi0 = np.array((1,)*Ndim)
        true_vac = np.array(optimize.minimize(v_, phi0, bounds = bnd).x)
        false_vac = [0]*Ndim
        points.append(Point(potential, delta, true_vac, false_vac))
        delta += step

    return points

def main():
    test_tag = sys.argv[1]
    Ndim = int(sys.argv[2])
    delta_min = float(sys.argv[3])
    delta_max = float(sys.argv[4])
    step = float(sys.argv[5])

    points = gen_points(Ndim, delta_min, delta_max, step)

    for point in points:
        T = 0
        v_ = lambda phi: point.potential.Vtot(phi,T)
        dv_ = lambda phi: point.potential.gradV(phi,T)
        path_pts = [point.true_vac, point.false_vac]

        blockPrint()
        tstart = time.process_time()
        instanton = pathDeformation.fullTunneling(path_pts, v_, dv_)
        tend = time.process_time()
        enablePrint()

        print('"{}",{},{},{}'.format(
            test_tag, point.delta, instanton.action, tend - tstart
        ))

    # print(gen_points(2, 0.1, 0.4, 0.01))

    # potential = PolynomialPotential(2, 0.4)
    
    # T = 0
    # v_ = lambda phi: potential.Vtot(phi,T)
    # dv_ = lambda phi: potential.gradV(phi,T)
    
    # bnd = ((.5, 3), (.5, 3))
    # phi0 = np.array((1,1))
    
    # true_vac = np.array(optimize.minimize(v_, phi0, bounds = bnd).x)
    # false_vac = np.array([0,0])
    # path_pts = [true_vac, false_vac]
    # print(path_pts)
    # instanton = pathDeformation.fullTunneling(path_pts, v_, dv_)
    # print(instanton.action)


if __name__ == "__main__":
    main()