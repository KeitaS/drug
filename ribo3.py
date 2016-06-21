#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True

@reaction_rules
def r30_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a > ~a | a * Lambda
    a + r30_u > r30_b | K_on * a * (r_u - r_min)
    r30_b > a + r30_u | K_off * r30_b

@reaction_rules
def r50_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a > ~a | a * Lambda
    a + r50_u > r50_b | K_on * a * (r_u - r_min)
    r50_b > a + r50_u | K_off * r50_b

@reaction_rules
def ribo_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a > ~a | a * Lambda
    a + r_u > r_b | K_on * a * (r_u - r_min)
    r_b > a + r_u | K_off * r_b
    r_b > a + r30_u + r50_u | Kd * r_u # r_bの解離


def createModel(r_max=65.8, r_min=19.3, K_D=1.0, K_t=6.1*10**-2, K_on=3.0, Lambda_0=1.35, Lambda_0_a=0.31, IC50=0.41, IC50_a=0.189, Kd=100., p=1., target=[]):
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
        if target[0] == "30s":
            r30_binding_reaction(a1_ex, a1, P_in, P_out, K_on, K_off, Lambda)
        elif target[0] == "50s":
            r50_binding_reaction(a1_ex, a1, P_in, P_out, K_on, K_off, Lambda)
        elif target[0] == "ribo":
            ribo_binding_reaction(a1_ex, a1, P_in, P_out, K_on, K_off, Lambda)
        
        if target[1] == "30s":
            r30_binding_reaction(a2_ex, a2, P_in, P_out, K_on, K_off, Lambda)
        elif target[1] == "50s":
            r50_binding_reaction(a2_ex, a2, P_in, P_out, K_on, K_off, Lambda)
        elif target[1] == "ribo":
            ribo_binding_reaction(a2_ex, a2, P_in, P_out, K_on, K_off, Lambda)

        ## ribo and subunit
        # production
        ~r30_u > r30_u | SUP
        ~r50_u > r50_u | SUP

        #bonding
        r30_u + r50_u > r_u | Ka * r30_u * r50_u

        # dissociation
        r_u > r30_u + r50_u | Kd * r_u
        r_u > ~r_u | r_u * Lambda
        r30_u > ~r30_u | r30_u * Lambda 
        r50_u > ~r50_u | r50_u * Lambda

        # diffusion
        r30_b > ~r30_b | r30_b * Lambda
        r50_b > ~r50_b | r50_b * Lambda
        r_b > ~r_b | r_b * Lambda

    return get_model()


def checkRiboDrug(drugNames=[]):
    """
    薬剤の標的を返す関数
    """
    targets = []
    for drugName in drugNames:
        if drugName == "30s":
            targets.append("30s")
        elif drugName == "50s":
            targets.append("50s")
        elif drugName == "ribo":
            targets.append("ribo")

    return targets


def run(drugs=[], step=10., legend=[], inpData={}, y0={"r30_u": 30., "r50_u": 30., "r_u": 30., "r_b": .0}):
    """
    ribosomeモデル実行関数

    drugs: [["drug_name": dose]]
    step: 
    legend: 
    inpData: 
    y0:

    ToDo:
    今後、drugの種類を判定する関数をつくり、この関数に導入する。
    """

    dataset = {"Lambda_0": 1.35, "Lambda_0_a": 0.31, "IC50": 0.41, "IC50_a": 0.189, "K_t": 6.1 * 10 ** -2, "r_min": 19.3}
    dNames = [] # drugの名前
    dDose = [] # drugのdose
    dTargets = [] # drugのターゲット
    
    dataset.update(inpData) # inpDataの内容をdatasetに入れる
    if drugs:
        for drug in drugs:
            dNames.append(drug[0])
            dDose.append(drug[1])
        dTargets = checkRiboDrug(dNames)
        dataset["target"] = dTargets # 標的リストをdatasetに入れる

    model = createModel(**dataset) # モデルを作成
    
    # y0に薬剤のDoseを入れる
    for index, dClass in enumerate(dNames):
        y0["a%d_ex" % (index + 1)] = dDose[index] # y0にdrugのdoseを入れる

    # legend
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
    p = 0.1

    ## drug data
    drug = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    Lambda_0 =  [1.35, 0.85, 0.40] # 1/h (1.35, 0.85, 0.40)
    Lambda_0_a = {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83} # 1/h 
    IC50 = {"Streptmycin": [0.41, 0.28, 0.196], "Kanamycin": [0.246, 0.096, 0.065], "Tetracycline": [0.5, 0.6, 1.45], "Chloramphenicol": [2.85, 2.65, 5.7]} # µg/ml 
    IC50_a = {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49} # µg/ml
    A_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20} # 

    xlabel = "Extracellular antibiotic concentration $a_{ex}$ ($\mu$M)"
    ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"


    # Kanamycinで比較
    # dose等はすべてKanamycinの値を用い、 
    name = "Kanamycin"
    dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "Lambda_0_a": Lambda_0_a[name], "IC50": IC50[name][0], "IC50_a": IC50_a[name], "p": p}
    legend = ["r_u"]
    
    drugs = [["30s", IC50[name][0]], ["30s", IC50[name][0]]]

    result, legend = run(drugs, inpData=dataset, legend=legend)
    savename = "20160620/%s_%s.png" % (drugs[0][0], drugs[1][0])
    makeGraph(np.array(result), savename, legend, name, xlabel, ylabel)

