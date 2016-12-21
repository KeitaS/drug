# coding: utf-8

import seaborn as sns
import matplotlib.pylab as plt
from ecell4 import *


# 定数
K_t = 6.1*10**-2

# 薬剤関連
K_D = 1. # 薬剤の解離定数と結合定数の比
K_on = 3.0 # 薬剤の結合定数

# リボソーム，growth関連
Lambda_0 = 1.35 # 初期培地でのGrowth Rate
r_min = 19.3 # リボソームの最小値
r_max = 65.8 # リボソームの最大値
r_u_0 = Lambda_0 / K_t + r_min # 定常状態でのr_uの値の予測

# リボソームサブユニット関連
Kd = 1. # リボソームサブユニットの解離定数
p = 1. # リボソームサブユニットの存在する比率
Ka = (Kd / K_t + r_u_0) * Lambda_0 / ((p * r_u_0) ** 2) # リボソームサブユニットの結合定数


with reaction_rules():
    r30_u + a1 == r30_b | (K_on, K_on)
    r50_u + a2 == r50_b | (K_on, K_on)
    r_u + a1 > r30_b + r50_u | K_on * (r_u - r_min)
    r_u + a2 > r50_b + r30_u | K_on * (r_u - r_min)
    r30_u + r50_u == r_u | (Ka, Kd * (r_u - r_min))
m = get_model()


def func(a1, a2, m):
    y0 = {"r_u": 30., "r30_u": 30., "r50_u": 30., "a1": a1, "a2": a2}
    return (run_simulation(10, model=m, y0=y0, species_list=["r_u"], return_type="array")[-1][1] - r_min) / r_max


plt.plot(range(0, 100), [func(d, 0, m) for d in range(0, 100)], "r-")
plt.plot(range(0, 100), [func(d * 0.5, d * 0.5, m) for d in range(0, 100)], "b-")
plt.show()


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
# plt.plot(range(0, 100), [func(d, 0) for d in range(0, 100)], "r-")
# plt.plot(range(0, 100), [func(d * 0.5, d * 0.5) for d in range(0, 100)], "b-")
# plt.show()
