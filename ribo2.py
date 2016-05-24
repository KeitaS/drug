#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True

"""
リボソームの動態
r50, r30, r
"""

def createModel(r_max=65.8, r_min=19.3, K_D=1.0, K_t=6.1*10**-2, K_on=3.0, Lambda_0=1.35, Lambda_0_a=0.31, IC50=0.41, IC50_a=0.189):
    Delta_r = r_max - r_min # µM
    K_off = K_on * K_D
    P_in = Delta_r * Lambda_0_a / 2.0 / IC50_a 
    P_out = (Lambda_0_a / 2) ** 2.0 / K_t / K_D

    with reaction_rules():
        Lambda = (r_u - r_min) * K_t
        Ka = Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))
        Kd = r_u * Lambda

        ~a_ex > a | P_in * a_ex
        a > ~a_ex | P_out * a
        a + r_u > r_b | K_on * a * (r_u - r_min)
        r_b > a + r_u | K_off * r_b
        # ~r_u > r_u | Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))
        # r_u > ~r_u | r_u * Lambda
        a > ~a | a * Lambda
        r_b > ~r_b | r_b * Lambda
        r30_u + r50_u == r_u | (Ka, Kd) # riboの合成

    return get_model()

def run(a_ex, step=10., legend=[], inpData={}, y0={"r30_u":60., "r50_u":60., "r_u": 30.0, "a": .0, "r_b": .0}):
    dataset = {"Lambda_0": 1.35, "Lambda_0_a": 0.31, "IC50": 0.41, "IC50_a": 0.189, "K_t": 6.1 * 10 ** -2, "r_min": 19.3}
    
    dataset.update(inpData)
    model = createModel(**dataset)
    y0["a_ex"] = a_ex

    if not legend:
        legend = y0.keys()

    runsim = run_simulation(step, solver="ode", y0=y0,
                            return_type="observer", model=model,
                            species_list=legend)
    data = runsim.data()
    return data, legend

def makeGraph(data, legend=[]):
    for i in range(len(data[0]) - 1):
        if legend[i]:
            plt.plot(data.T[0], data.T[i+1], label=legend[i])
        else:
            plt.plot(data.T[0], data.T[i+1])
    plt.legend(loc="upper right")
    plt.show()
    # plt.savefig("result/ribo2_2_3.png", dpi=200)

if __name__ == "__main__":
    dataset = {"Lambda_0": 0.982371812727}
    result, legend = run(.5, inpData=dataset)
    makeGraph(np.array(result), legend)
