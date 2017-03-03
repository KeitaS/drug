# coding: utf-8
"""
    r30, r50に，それぞれ2つ薬剤が結合できるようにしたもの．
    State : add_reaction
"""

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True
import itertools

# main
@reaction_rules
def r30_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda, r_min, num):
    ~a_ex == a | (P_in * a_ex, P_out)
    if num == 0: # 一つ目の薬剤の場合
        a + r30_u_u == r30_b_u | (K_on, K_off)
        a + r30_u_b == r30_b_b | (K_on, K_off)
        a + r_u > r30_b_u + r50_u_u | K_on * a * (r_u - r_min) # add
        r30_b_u > ~r30_b_u | Lambda * r30_b_u
    elif num == 1: # 二つ目の薬剤の場合
        a + r30_u_u == r30_u_b | (K_on, K_off)
        a + r30_b_u == r30_b_b | (K_on, K_off)
        a + r_u > r30_u_b + r50_u_u | K_on * a * (r_u - r_min) # add
        r30_u_b > ~r30_u_b | Lambda * r30_u_b
    a > ~a | Lambda * a

@reaction_rules
def r50_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda, r_min, num):
    ~a_ex == a | (P_in * a_ex, P_out)
    if num == 0: # 一つ目の薬剤の場合
        a + r50_u_u == r50_b_u | (K_on, K_off)
        a + r50_u_b == r50_b_b | (K_on, K_off)
        a + r_u > r50_b_u + r30_u_u | K_on * a * (r_u - r_min) # add
        r50_b_u > ~r50_b_u | Lambda * r50_b_u
    elif num == 1: # 二つ目の薬剤の場合
        a + r50_u_u == r50_u_b | (K_on, K_off)
        a + r50_b_u == r50_b_b | (K_on, K_off)
        a + r_u > r50_u_b + r30_u_u | K_on * a * (r_u - r_min) # add
        r50_u_b > ~r50_u_b | Lambda * r50_u_b
    a > ~a | Lambda * a

@reaction_rules
def ribo_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda, r_min, num):
    ~a_ex == a | (P_in * a_ex, P_out)
    a + r_u == r_b | (K_on * a * (r_u - r_min), K_off)
    r_b > a + r30_u_u + r50_u_u | Kd
    a > ~a | Lambda * a

def createModel(drugs=[], dataset={}):
    # 定数
    K_t = 6.1 * 10 ** -2

    # 薬剤関連
    K_D = 1. # 薬剤の解離定数と結合定数の比
    K_on = 3.0 # 薬剤の結合定数
    K_off = K_on * K_D # 薬剤の解離定数

    # リボソーム，growth関連
    Lambda_0 = 1.35 # 初期培地でのGrowth Rate
    r_min = 19.3 # リボソームの最小値
    r_max = 65.8 # リボソームの最大値
    Delta_r = r_max - r_min # リボソームの最大値と最小値の差
    r_u_0 = Lambda_0 / K_t + r_min # 定常状態でのr_uの値の予測

    # リボソームサブユニット関連
    Kd = 1. # リボソームサブユニットの解離定数
    p = 1. # リボソームサブユニットの存在する比率
    Ka = (Kd / K_t + r_u_0) * Lambda_0 / ((p * r_u_0) ** 2) # リボソームサブユニットの結合定数

    # 引数でdatasetを受け取った場合，更新する
    if dataset:
        for key, value in dataset.items():
            exec("{} = {}".format(key, value))

    with reaction_rules():
        # expression
        Lambda = (r_u - r_min) * K_t
        SUP = (Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1 / K_t / Delta_r))) * (1 + p) # subunit product expression

        # drug reaction
        if drugs:
            for index, drug in enumerate(drugs):
                a_ex = _eval("a{}_ex".format(index))
                a = _eval("a{}".format(index))
                P_in = (r_max - r_min) * drug["Lambda_0_a"] / 2.0 / drug["IC50_a"] # 薬剤の流入
                P_out = (drug["Lambda_0_a"] / 2) ** 2.0 / K_t / K_D # 薬剤の流出
                exec('{}_binding_reaction(a_ex, a, P_in, P_out,\
                                           K_on, K_off, Lambda, r_min, index)'\
                                           .format(drug["target"]))

        # another reaction
        ~r30_u_u > r30_u_u | SUP
        ~r50_u_u > r50_u_u | SUP
        r30_u_u + r50_u_u == r_u | (Ka, Kd * (r_u - r_min))

        r30_u_u > ~r30_u_u | Lambda * r30_u_u
        r50_u_u > ~r50_u_u | Lambda * r50_u_u
        r_u > ~r_u | Lambda * r_u

        r30_b_b > ~r30_b_b | Lambda * r30_b_b
        r50_b_b > ~r50_b_b | Lambda * r50_b_b
        r_b > ~r_b | Lambda * r_b

    return get_model()

def run(drugs=[], step=50., legend=[], inpData={}, y0={"r30_u_u": 30., "r50_u_u": 30., "r_u": 30.}):
    model = createModel(drugs, inpData)
    #
    # for i, rr in enumerate(model.reaction_rules()):
    #     print("{}, {}".format(i, rr.as_string()))

    if drugs:
        for index, drug in enumerate(drugs):
            y0["a{}_ex".format(index)] = drug["dose"]

    if not legend:
        legend = list(y0.keys())
    runsim = run_simulation(step, solver="ode", y0=y0,
                            return_type="observer", model=model,
                            species_list=legend)
    data = runsim.data()
    # legend.insert(0, "t")
    # data = pd.DataFrame(runsim.data(), columns=legend)
    return data, legend


# sub modules
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
    types = {"Streptmycin": "r30", "Kanamycin": "r30", "Tetracycline": "r30", "Chloramphenicol": "r50"}

    drugData = {"name": drugName,
                "target": types[drugName],
                "dose": .0,
                "Lambda_0_a": Lambda_0_a[drugName],
                "IC50": IC50[drugName][medium],
                "IC50_a": IC50_a[drugName],
                "K_ma": 15.
                }

    return drugData

def calcGrowthRate(r_u, r_min=19.3, K_t=6.1*10**-2, Lambda_0=1.35):
    result = (r_u - r_min) * K_t / Lambda_0
    return result

def calcIC(dNames, a_ex, target):
    """
    二分法でIC〜〜を計算
    """
    calc_result = {}
    for dName in dNames:
        print(dName)
        drugs = [makeDrugDatas(dName, 0)] # 薬剤データの作成
        dose_max = a_ex[dName] * 3
        dose_min = .0
        result = 1.
        dose = .0
        count = 0
        while abs(target - abs(result)) > 0.01 and count <= 10:
            before_dose = dose
            dose = (dose_max + dose_min) / 2.
            if dose == before_dose:
                count += 1
            print(dose)
            print(count)
            drugs[0]["dose"] = dose
            result, legend = run(drugs, legend=["r_u"])
            result = calcGrowthRate(result[-1][1])
            if result < target:
                dose_max = dose
            else:
                dose_min = dose
        calc_result[dName] = dose
    return calc_result
