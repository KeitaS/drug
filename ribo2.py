#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True

"""
リボソームの動態
r50, r30, r
"""

def createModel():
    Ka = 0.01
    Kd = 0.3
    with reaction_rules():
        r30_u + r50_u == r_u | (Ka, Kd) # riboの合成

    return get_model()

def run(step=10., legend=[], dataset={}):
    model = createModel()
    y0 = {"r30_u":60., "r50_u":60., "r_u": 0.}
    if dataset:
        y0.update(dataset)

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
    # plt.savefig("result/test.png", dpi=200)
    plt.legend(loc="upper right")
    plt.show()

if __name__ == "__main__":
    result, legend = run()
    makeGraph(np.array(result), legend)
