#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True


def createModel(r_max=65.8, r_min=19.3, K_D=1.0, K_t=6.1*10**-2, K_on=3.0, Lambda_0=1.35, Lambda_0_a=0.31, IC50=0.41, IC50_a=0.189, Kd=100., p=1.):
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
    
    Consideration:
    薬剤固有の値はIC50_a, IC50, Lambda_0_a
    複数入ってきた時の対応はどうするか
    """

    Delta_r = r_max - r_min # µM
    K_off = K_on * K_D # riboと薬剤との結合
    P_in = Delta_r * Lambda_0_a / 2.0 / IC50_a # 薬剤の流入
    P_out = (Lambda_0_a / 2) ** 2.0 / K_t / K_D # 薬剤の流出
    r_u_0 = Lambda_0 / K_t + r_min
    Ka = (Lambda_0 + Kd) / (p * p * r_u_0)

    with reaction_rules():
        ### expression
        Lambda = (r_u - r_min) * K_t
        SUP = (Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))) * (1 + p) # subunit production expression
        

        ### reaction
        ## drug
        # a1(r30)
        ~a1_ex > a1 | P_in * a1_ex
        a1 > ~a1_ex | P_out * a1
        a1 > ~a1 | a1 * Lambda
        
        # a2(r50)
        ~a2_ex > a2 | P_in * a2_ex
        a2 > ~a2_ex | P_out * a2
        a2 > ~a2 | a2 * Lambda

        # a3(r)
        ~a3_ex > a3 | P_in * a3_ex
        a3 > ~a3_ex | P_out * a3
        a3 > ~a3 | a3 * Lambda

        ## ribo and subunit
        r30_u + r50_u > r_u | Ka * r30_u * r50_u # bonding
        r_u > r30_u + r50_u | Kd * r_u # dissociation
        ~r30_u > r30_u | SUP # production
        ~r50_u > r50_u | SUP # production
        r_u > ~r_u | r_u * Lambda # diffusion
        r30_u > ~r30_u | r30_u * Lambda # diffusion
        r50_u > ~r50_u | r50_u * Lambda # diffusion

        ## bond drug
        # r30
        a1 + r30_u > r30_b | K_on * a1 * (r_u - r_min)
        r30_b > a1 + r30_u | K_off * r30_b
        r30_b > ~r30_b | r30_b * Lambda

        # r50
        a2 + r50_u > r50_b | K_on * a2 * (r50_u - r_min)
        r50_b > a2 + r50_u | K_off * r50_b
        r50_b > ~r50_b | r50_b * Lambda
        
        # r
        a3 + r_u > r_b | K_on * a3 * (r_u - r_min)
        r_b > a3 + r_u | K_off * r_b
        r_b > a3 + r30_u + r50_u | Kd * r_u # r_bの解離
        r_b > ~r_b | r_b * Lambda

    return get_model()

def run(a1_ex=.0, a2_ex=.0, a3_ex=.0, step=10., legend=[], inpData={}, y0={"r30_u":30., "r50_u":30., "r_u":30., "a": .0, "r_b": .0}):
    dataset = {"Lambda_0": 1.35, "Lambda_0_a": 0.31, "IC50": 0.41, "IC50_a": 0.189, "K_t": 6.1 * 10 ** -2, "r_min": 19.3}
    
    dataset.update(inpData)
    model = createModel(**dataset)
    y0["a1_ex"] = a1_ex # target S30
    y0["a2_ex"] = a2_ex # target S50
    y0["a3_ex"] = a3_ex # target ribo

    if not legend: 
        legend = y0.keys()

    runsim = run_simulation(step, solver="ode", y0=y0,
                            return_type="observer", model=model,
                            species_list=legend)
    data = runsim.data()
    return data, legend


def makeGraph(data, savename, legend=[], title="", xlabel="", ylabel=""):
    """
    グラフ作成用モジュール
    """
    for i in range(len(data[0]) - 1):
        if legend[i]:
            plt.plot(data.T[0], data.T[i+1], label=legend[i])
        else:
            plt.plot(data.T[0], data.T[i+1])

    if title: # titleがあったらつける
        plt.title(title)

    if xlabel: # xlabelがあったらつける
        plt.xlabel(xlabel)

    if ylabel: # ylabelがあったらつける
        plt.ylabel(ylabel)

    plt.legend(loc="upper right")
    # plt.show()
    plt.savefig("result/%s" % (savename), dpi=200)
    plt.close()


if __name__ == "__main__":
    r_min = 19.3
    K_t = 6.1 * 10 ** -2
    Kd = 5.
    p = 0.05

    ## drug data
    drug = ["streptmycin", "kanamycin", "tetracycline", "chloramphenicol"]
    Lambda_0 =  [1.35, 0.85, 0.40] # 1/h (1.35, 0.85, 0.40)
    Lambda_0_a = {"streptmycin": 0.31, "kanamycin": 0.169, "tetracycline": 5.24, "chloramphenicol": 1.83} # 1/h 
    IC50 = {"streptmycin": [0.41, 0.28, 0.196], "kanamycin": [0.246, 0.096, 0.065], "tetracycline": [0.5, 0.6, 1.45], "chloramphenicol": [2.85, 2.65, 5.7]} # µg/ml 
    IC50_a = {"streptmycin": 0.189, "kanamycin": 0.05, "tetracycline": 0.229, "chloramphenicol": 2.49} # µg/ml
    a_ex = {"streptmycin": 0.6, "kanamycin": 0.5, "tetracycline":2, "chloramphenicol": 20} # 

    xlabel = "Extracellular antibiotic concentration $a_{ex}$ ($\mu$M)"
    ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"


    # カナマイシンで比較
    name = "kanamycin"
    dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "Lambda_0_a": Lambda_0_a[name], "IC50": IC50[name][0], "IC50_a": IC50_a[name], "p": p}
    legend = ["r_u"]
    data = []

    for i in np.linspace(0, a_ex[name], 51):
        result, legend = run(a3_ex=i, inpData=dataset, legend=legend)
        data.append([i, (result[-1][1] - r_min) * K_t / Lambda_0[0]])

    savename = "20160614/ribo2/%s_2.png" % (name) 
    makeGraph(np.array(data), savename, legend, name, xlabel, ylabel)

