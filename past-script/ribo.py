#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import pandas as pd
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True


def createModel(r_max=65.8, r_min=19.3, K_D=0.1, K_t=6.1*10**-2, K_on=60.0, Lambda_0=1.35, Lambda_0_a=0.31, IC50=0.41, IC50_a=0.189):
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
        Lambda = (r_u - r_min) * K_t

        ~a_ex > a | P_in * a_ex
        a > ~a_ex | P_out * a
        a + r_u > r_b | K_on * a * (r_u - r_min)
        r_b > a + r_u | K_off * r_b
        ~r_u > r_u | Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))
        a > ~a | a * Lambda
        r_u > ~r_u | r_u * Lambda
        r_b > ~r_b | r_b * Lambda

    return get_model()

def sim(drug, a_ex, inpData={"Lambda_0": 1.35}):
    r_min = 19.3
    K_t = 6.1*10**-2
    Lambda_0 = inpData["Lambda_0"]
    data = [];
    for dose in a_ex:
        drug["dose"] = dose
        result, legend = run(drug, inpData=inpData)
        data.append((result[-1][1] - r_min) * K_t / Lambda_0)
    return data

def run(drug, step=50., legend=["r_u"], inpData={}, y0={"a": .0, "r_u": 30.0, "r_b": .0}):
    dataset = {"Lambda_0": 1.35, "Lambda_0_a": drug["Lambda_0_a"], "IC50": drug["IC50"], "IC50_a": drug["IC50_a"], "K_t": 6.1 * 10 ** -2, "r_min": 19.3}

    dataset.update(inpData)
    model = createModel(**dataset)
    y0["a_ex"] = drug["dose"]

    if not legend:
        legend = y0.keys()
    runsim = run_simulation(step, solver="ode", y0=y0,
                            return_type="observer", model=model,
                            species_list=legend)
    data = runsim.data()
    return data, legend

def makeDrugDatas(drugName, medium=0):
    """
    薬剤データを作成して返す関数
    drugName : 投与する薬剤の名前
    medium : 培地の番号(0: ,
                      1: ,
                      2: ,
                     )
    """
    Lambda_0_a = {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83} # 1/h
    IC50 = {"Streptmycin": [0.41, 0.28, 0.196], "Kanamycin": [0.246, 0.096, 0.065], "Tetracycline": [0.5, 0.6, 1.45], "Chloramphenicol": [2.85, 2.65, 5.7]} # µg/ml
    IC50_a = {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49} # µg/ml
    # types = {"Streptmycin": "r30", "Kanamycin": "r30", "Tetracycline": "r30", "Chloramphenicol": "r50"}

    drugData = {"name": drugName,
                # "target": types[drugName],
                "dose": .0,
                "Lambda_0_a": Lambda_0_a[drugName],
                "IC50": IC50[drugName][medium],
                "IC50_a": IC50_a[drugName]
                }

    return drugData

if __name__ == "__main__":
    r_min = 19.3
    K_t = 6.1 * 10 ** -2

    ## drug data
    drugName = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    Lambda_0 =  [1.35, 0.85, 0.40] # 1/h (1.35, 0.85, 0.40)
    A_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}

    xlabel = "Extracellular antibiotic concentration $a_{ex}$ ($\mu$M)"
    ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"
    legend = ["$Gly$", "$Gly_{CAA}$", "$Gly_{RDM}$"]
    color = ["g", "b", "r"]

    doses = np.linspace( 0, A_ex[drugName[0]], 100)
    drug = makeDrugDatas(drugName[0])
    result = sim(drug, doses)
    plt.plot(doses, result)
    plt.show()
