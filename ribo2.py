#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True

"""
リボソームの動態
r50, r30, r
"""

def createModel(r_max=65.8, r_min=19.3, K_D=1.0, K_t=6.1*10**-2, K_on=3.0, Lambda_0=1.35, Lambda_0_a=0.31, IC50=0.41, IC50_a=0.189, Ka=3., Kd=.0):
    """
    リボソームモデルを構成するモジュール
    r_max: µM
    r_min: µM
    K_D: none
    K_t: 1/µM/h
    K_on: ???
    K_off: ???
    Delta_r: µM
    Lambda_0: default medium = Gly_RDM
    Lambda_0_a: default drug = streptmycin
    IC50: default medium = Gly_RDM, default drug = streptmycin
    IC50_a: default drug = streptmycin
    Ka: subunitの結合
    Kd: subunitの解離

    薬剤によって変更する必要がある値:
    Lambda_0_a, IC50, IC50_a

    培地によって変更する値(COBRAの結果):
    Lambda_0
    """

    Delta_r = r_max - r_min # µM
    K_off = K_on * K_D
    P_in = Delta_r * Lambda_0_a / 2.0 / IC50_a 
    P_out = (Lambda_0_a / 2) ** 2.0 / K_t / K_D

    with reaction_rules():
        ## expression
        Lambda = (r_u - r_min) * K_t
        SUP = Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r)) # subunit production expression

        ## reaction
        # a
        ~a_ex > a | P_in * a_ex
        a > ~a_ex | P_out * a
        a > ~a | a * Lambda
        
        # ribo and subunit
        r30_u + r50_u > r_u | Ka * r30_u * r50_u
        r_u > r30_u + r50_u | Kd * r_u 
        ~r30_u > r30_u | sup # production
        ~r50_u > r50_u | sup # production
        r_u > ~r_u | r_u * Lambda # diffusion
        r30_u > ~r30_u | r30_u * Lambda # diffusion
        r50_u > ~r50_u | r50_u * Lambda # diffusion

        # bond drug
        a + r_u > r_b | K_on * a * (r_u - r_min)
        r_b > a + r_u | K_off * r_b
        r_b > ~r_b | r_b * Lambda

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

def makeGraph(data, savename, legend=[]):
    for i in range(len(data[0]) - 1):
        if legend[i]:
            plt.plot(data.T[0], data.T[i+1], label=legend[i])
        else:
            plt.plot(data.T[0], data.T[i+1])
    plt.legend(loc="upper right")
    # plt.show()
    plt.savefig("result/%s" % (savename), dpi=200)
    plt.close()

if __name__ == "__main__":
    """
    savename = "ribo2_1.png"
    dataset = {"Lambda_0": 0.982371812727}
    legend = ["r_u", "r_b", "a_ex", "a"]
    result, legend = run(.0, inpData=dataset, legend=legend)
    makeGraph(np.array(result), savename, legend)
    """
    count = 0
    for i in np.linspace(0, 1, 101):
        legend = ["r_u", "r_b", "a_ex", "a"]
        dataset = {"Ka": i, "Lambda_0": 0.982371812727}
        savename = "ribo2_%d.png" % (count)
        result, legend = run(.0, inpData=dataset, legend=legend)
        makeGraph(np.array(result), savename, legend)
        print "Ka: %f finish." % (i)
        count += 1
