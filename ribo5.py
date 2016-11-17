#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True

@reaction_rules
def r30_binding_reaction(a_ex, a, r30_b, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a + r30_u > r30_b | K_on * a * r30_u
    r30_b > a + r30_u | K_off * r30_b
    a > ~a | a * Lambda # dilution
    r30_b > ~r30_b | r30_b * Lambda # dilution


@reaction_rules
def r50_binding_reaction(a_ex, a, r50_b, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a + r50_u > r50_b | K_on * a * r50_u
    r50_b > a + r50_u | K_off * r50_b # dissociation
    a > ~a | a * Lambda # dilution
    r50_b > ~r50_b | r50_b * Lambda # dilution


@reaction_rules
def ribo_binding_reaction(a_ex, a, r_b, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a + r_u > r_b | K_on * a * (r_u - r_min)
    r_b > a + r_u | K_off * r_b
    r_b > a + r30_u + r50_u | Kd * r_b # dissociation
    a > ~a | a * Lambda # dilution
    r_b > ~r_b | r_b * Lambda # dilution


def createModel(drugs=[], r_max=65.8, r_min=19.3, K_D=1., K_t=6.1*10**-2, K_on=3.0, Lambda_0=1.35, Kd=1., p=1., K_ma=3.):
    """
    リボソームモデルを構成するモジュール
    r_max: µM
    r_min: µM
    K_D: none : 薬剤の結合定数
    K_t: 1/µM/h
    K_on: ??? : 結合定数
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
    r_uは、r_min以下の値の時は機能していないことになっている。
    そのため、r_uの機能する部分がr30_uとr50_uに分離するようにモデルを改変した。
    r_uの考慮により、Ka、SUPも付随して変更が起こっている。

    式:
        r_u : r30_u = 1 : p
        Ka * r30_u * r50_u - Kd * (r_u - r_min) = r_u * Lambda
        r_u - r_min = Lambda / K_t
    """

    Delta_r = r_max - r_min # µM
    K_off = K_on * K_D # riboと薬剤との結合
    r_u_0 = Lambda_0 / K_t + r_min # 定常の時にこの値になる。
    Ka = (Kd / K_t + r_u_0) * Lambda_0 / ((p * r_u_0) ** 2)
    K_ma1 = K_ma
    K_ma2 = K_ma

    with reaction_rules():
        ### expression
        Lambda = (r_u - r_min) * K_t
        SUP = (Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))) * (1 + p) # subunit product expression

        ### reaction
        ## drug
        if len(drugs) > 0:
            K_on_exp = K_on * 1 / (1 + a2 / K_ma2)
            # K_off_exp = K_on_exp * K_D
            K_off_exp = K_off
            P_in = drugs[0]["P_in"] * 1 / (1 + a2 / K_ma2)

            if drugs[0]["type"] == "30s":
                # r30_binding_reaction(a1_ex, a1, r30_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
                # r30_binding_reaction(a1_ex, a1, r30_1_b, P_in, drugs[0]["P_out"], K_on, K_off, Lambda)
                r30_binding_reaction(a1_ex, a1, r30_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on_exp, K_off_exp, Lambda)
            elif drugs[0]["type"] == "50s":
                # r50_binding_reaction(a1_ex, a1, r50_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
                # r50_binding_reaction(a1_ex, a1, r50_1_b, P_in, drugs[0]["P_out"], K_on, K_off, Lambda)
                r50_binding_reaction(a1_ex, a1, r50_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on_exp, K_off_exp, Lambda)
            elif drugs[0]["type"] == "ribo":
                # ribo_binding_reaction(a1_ex, a1, r_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
                # ribo_binding_reaction(a1_ex, a1, r_1_b, P_in, drugs[0]["P_out"], K_on, K_off, Lambda)
                ribo_binding_reaction(a1_ex, a1, r_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on_exp, K_off_exp, Lambda)

        if len(drugs) > 1:
            K_on_exp = K_on * 1 / (1 + a1 / K_ma1)
            # K_off_exp = K_on_exp * K_D
            K_off_exp = K_off
            P_in = drugs[1]["P_in"] * 1 / (1 + a1 / K_ma1)

            if drugs[1]["type"] == "30s":
                # r30_binding_reaction(a2_ex, a2, r30_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on, K_off, Lambda)
                # r30_binding_reaction(a2_ex, a2, r30_2_b, P_in, drugs[1]["P_out"], K_on, K_off, Lambda)
                r30_binding_reaction(a2_ex, a2, r30_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on_exp, K_off_exp, Lambda)
            elif drugs[1]["type"] == "50s":
                # r50_binding_reaction(a2_ex, a2, r50_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on, K_off, Lambda)
                # r50_binding_reaction(a2_ex, a2, r50_2_b, P_in, drugs[1]["P_out"], K_on, K_off, Lambda)
                r50_binding_reaction(a2_ex, a2, r50_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on_exp, K_off_exp, Lambda)
            elif drugs[1]["type"] == "ribo":
                # ribo_binding_reaction(a2_ex, a2, r_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on, K_off, Lambda)
                # ribo_binding_reaction(a2_ex, a2, r_2_b, P_in, drugs[1]["P_out"], K_on, K_off, Lambda)
                ribo_binding_reaction(a2_ex, a2, r_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on_exp, K_off_exp, Lambda)

        ## ribo and subunit
        # production
        ~r30_u > r30_u | SUP
        ~r50_u > r50_u | SUP

        # bonding
        r30_u + r50_u > r_u | Ka * r30_u * r50_u

        # dissociation
        r_u > r30_u + r50_u | Kd * (r_u - r_min)

        # dilution
        r_u > ~r_u | r_u * Lambda
        r30_u > ~r30_u | r30_u * Lambda
        r50_u > ~r50_u | r50_u * Lambda

    return get_model()


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
    types = {"Streptmycin": "30s", "Kanamycin": "30s", "Tetracycline": "30s", "Chloramphenicol": "50s"}

    drugData = {"name": drugName, "type": types[drugName], "dose": .0, "Lambda_0_a": Lambda_0_a[drugName], "IC50": IC50[drugName][medium], "IC50_a": IC50_a[drugName]}

    return drugData


def run(drugs=[], step=50., legend=[], inpData={}, y0={"r30_u": 30., "r50_u": 30., "r_u": 30., "r_b": .0}):
    """
    ribosomeモデル実行関数

    drugs: [{"name":, "type":, "dose":, }]
    step:
    legend:
    inpData: createModelに渡すDataset
    y0:
    """

    dataset = {"r_max": 65.8, "r_min": 19.3, "K_D": 1.0, "K_t": 6.1*10**-2, "K_on": 3.0, "Lambda_0": 1.35, "Kd": 1., "p": 1., "K_ma": 3.} # 基本のデータセット
    dataset.update(inpData) # inpDataの内容をdatasetに入れる

    # P_in, P_outをdrugのデータに入れる
    for drug in drugs: # drug dataに基づく情報を用いた
        drug["P_in"] = (dataset["r_max"] - dataset["r_min"]) * drug["Lambda_0_a"] / 2.0 / drug["IC50_a"] # 薬剤の流入
        drug["P_out"] = (drug["Lambda_0_a"] / 2) ** 2.0 / dataset["K_t"] / dataset["K_D"] # 薬剤の流出

    if drugs:
        dataset["drugs"] = drugs
        # y0に薬剤のdoseを入れる
        for index, drug in enumerate(drugs):
            y0["a%d_ex" % (index + 1)] = drug["dose"]

    model = createModel(**dataset) # モデルを作成

    # legend
    if not legend:
        legend = y0.keys()

    runsim = run_simulation(step, solver="ode", y0=y0,
                            return_type="observer", model=model,
                            species_list=legend)
    data = runsim.data()
    return data, legend


def calcGrowthrate(a, r_min=19.3, K_t=6.1*10**-2, Lambda_0=1.35):
    """
    growth_rateを計算する関数
    """
    result = (a - r_min) * K_t / Lambda_0
    return result


def checkLinerType(inp, eps_rel, pattern, buffpoint=0.):
    """
    inp: growth rate list
    eps_rel: 差の許容率
    pattern: linertypeを計算するときのパターン
    """
    inp = np.array(inp)
    dinp = inp[1:] - inp[:-1] # 差のリスト
    before_type = 0 # -1: <, 0: ==, 1: >
    linertype = 0 # -1: synergistic, 0: additive, 2: antagonistic`

    if pattern == 0:
        for i in range(len(dinp)):
            if abs(dinp[i]) > inp[i] * eps_rel: #
                if dinp[i] < 0:
                    current_type = -1 # 減少
                else:
                    current_type = 1 # 増加
            else:
                current_type = 0

            if i == 0:
                before_type = current_type
            else:
                if current_type == 0 or before_type == 0 or current_type == before_type:
                    pass
                else:
                    linerList = np.linspace(inp[0], inp[-1], len(inp)) # 両端点を線形で結んだList
                    linertype = linerList[i] - inp[i] #

                    return linertype
                    break

    elif pattern == 1:
        upper_bound = max(inp[0], inp[-1])
        lower_bound = min(inp[0], inp[-1])
        max_inp = max(inp)
        min_inp = min(inp)
        if max_inp > upper_bound: # antagonistic
            linertype = upper_bound - max_inp
        elif min_inp < lower_bound: # synergistic
            linertype = lower_bound - min_inp

    elif pattern == 2:
        if buffpoint != 0:
            upper_bound = max(inp[0], inp[-1])
            lower_bound = min(inp[0], inp[-1])
            max_inp = max(inp)
            min_inp = min(inp)
            if max_inp > upper_bound: # antagonistic
                linertype = -(max_inp / buffpoint)
            elif min_inp < lower_bound: # synergistic
                linertype = lower_bound - min_inp

    return linertype


def calcIC(dNames, a_ex, target):
    """
    二分法でIC〜〜を計算
    """
    calc_result = {}
    for dName in dNames:
        drugs = [makeDrugDatas(dName, 0)] # 薬剤データの作成
        dose_max = a_ex[dName] * 3
        dose_min = .0
        result = 1.
        dose = .0
        while abs(target - abs(result)) > 0.01:
            dose = (dose_max + dose_min) / 2.
            drugs[0]["dose"] = dose
            result, legend = run(drugs, step=100, legend=["r_u"])
            result = calcGrowthrate(result[-1][1])
            if result < target:
                dose_max = dose
            else:
                dose_min = dose
        calc_result[dName] = dose
    return calc_result


def calcBufferingPoint(dNames, doses):
    """
    二分法でBufferingPointを計算
    """
    a_range = [0, doses[0]]
    a_mid = 0
    b_range = [0, doses[1]]
    b_mid = doses[1]
    result = 1
    drugA = [makeDrugDatas(dNames[0])]
    drugB = [makeDrugDatas(dNames[1])]
    result_a = 0
    result_b = 0
    while abs(result) > 0.01:
        if result > 0:
            a_range[0] = a_mid
            b_range[1] = b_mid
        else:
            a_range[1] = a_mid
            b_range[0] = b_mid
        a_mid = np.linspace(a_range[0], a_range[1], 3)[1]
        b_mid = np.linspace(b_range[0], b_range[1], 3)[::-1][1]
        drugA[0]["dose"] = a_mid
        drugB[0]["dose"] = b_mid
        result, legend = run(drugA, step=100, legend=["r_u"])
        result_a = calcGrowthrate(result[-1][1])
        result, legend = run(drugB, step=100, legend=["r_u"])
        result_b = calcGrowthrate(result[-1][1])
        result = result_a - result_b

    return min(result_a, result_b)


if __name__ == "__main__":
    # import
    import itertools as itr
    import seaborn as sns
    import pandas as pd
    from modules import *

    # 保存用ディレクトリの作成
    savedir = "./images/ribo5"
    makedir(savedir)

    ## drug data
    dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}

    # IC30 = calcIC(dNames, a_ex, .3) # IC30を計算
    IC30 = {'Kanamycin': 0.6761398315429688, 'Streptmycin': 1.4652465820312497, 'Chloramphenicol': 22.5, 'Tetracycline': 5.25}
