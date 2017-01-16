from ecell4 import *
import numpy as np
import matplotlib.pylab as plt


pin, pout = 1., 1.
ka, kd = 3., 3.
k1, k2 = 0.05, 1.0
rmin, rmax = 19.3, 65.8

with reaction_rules():
    R30 + a == R30_a | (ka * R30 * a, kd)
    R50 + b == R50_b | (ka * R50 * b, kd)
    R50 + R30 == R | (k1, k2 * (R - rmin))
    R + a > R30_a + R50 | ka * (R - rmin)
    R + b > R30 + R50_b | ka * (R - rmin)

m = get_model()

def run(m, y0):
    return run_simulation(30.0, model=m, y0=y0, species_list=['R'], return_type='array')[-1][1]

x = np.linspace(0., 100., 51)
# phi = run(m, {"R30": 20., "R50": 20.})
# plt.plot(x, [run(m, {"R30": 20., "R50": 20., "a": xi * pin / pout}) / phi for xi in x], "--", label="single")
# plt.plot(x, [run(m, {"R30": 20., "R50": 20., "a": xi * pin / pout * 0.5, "b": xi * pin / pout * 0.5}) / phi for xi in x], label="double")
phi = run(m, {"R": rmax})
plt.plot(x, [run(m, {"R": rmax, "a": xi * pin / pout}) / phi for xi in x], label="single")
plt.plot(x, [run(m, {"R": rmax, "a": xi * pin / pout * 0.5, "b": xi * pin / pout * 0.5}) / phi for xi in x], label="double")
plt.legend(loc="best")
plt.xlabel("dose")
plt.ylabel("Normalized Growth Rate $\lambda / \lambda_0$")
plt.savefig("images/simple/antagonistic.png", dpi=300)

###
# k1 = k2 = k3 = k4 = 3.0
# k5, k6 = 0.05, 1.0
# R_min = 19.3
# R_max = 65.8
#
# with reaction_rules():
#     R30 + a == R30_a | (k1, k2)
#     R50 + b == R50_b | (k3, k4)
#     R + a > R50 + R30_a | k1 * (R - R_min)
#     R + b > R30 + R50_b | k3 * (R - R_min)
#     R30 + R50 == R | (k5, k6 * (R - R_min))
#
# m = get_model()
#
# def func(a, b):
#     return (run_simulation(10, model=m, y0={"R": R_max, "a": a, "b": b}, species_list=["R"], return_type="array")[-1][1] - R_min) / R_max
#
# plt.plot(range(0, 100), [func(d, 0) for d in range(0, 100)], "r-", label="single")
# plt.plot(range(0, 100), [func(d * 0.5, d * 0.5) for d in range(0, 100)], "b-", label="double")
# plt.legend(loc="best")
# plt.show()
