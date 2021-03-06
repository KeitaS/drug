#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools as itr
import sys

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


def createModel(drugs=[], r_max=65.8, r_min=19.3, K_D=1., K_t=6.1*10**-2, K_on=3.0, Lambda_0=1.35, Kd=1., p=1., modif=0):
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

    with reaction_rules():
        ### expression
        Lambda = (r_u - r_min) * K_t
        SUP = (Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))) * (1 + p) # subunit product expression

        ### reaction
        ## drug
        if len(drugs) > 0:
            # "単剤と多剤で分ける
            if len(drugs) == 1:
                K_on_exp = K_on
                P_in = drugs[0]["P_in"]
            else:
                K_on_exp = K_on * 1 / (1 + a2 / drugs[1]["K_ma"])
                P_in = drugs[0]["P_in"] * 1 / (1 + a2 / drugs[1]["K_ma"])
            # K_off_exp = K_on_exp * K_D
            K_off_exp = K_off


            if drugs[0]["type"] == "30s":
                if modif == 0: r30_binding_reaction(a1_ex, a1, r30_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
                elif modif == 1: r30_binding_reaction(a1_ex, a1, r30_1_b, P_in, drugs[0]["P_out"], K_on, K_off, Lambda)
                elif modif == 2:r30_binding_reaction(a1_ex, a1, r30_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on_exp, K_off_exp, Lambda)
            elif drugs[0]["type"] == "50s":
                if modif == 0: r50_binding_reaction(a1_ex, a1, r50_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
                elif modif == 1: r50_binding_reaction(a1_ex, a1, r50_1_b, P_in, drugs[0]["P_out"], K_on, K_off, Lambda)
                elif modif == 2:r50_binding_reaction(a1_ex, a1, r50_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on_exp, K_off_exp, Lambda)
            elif drugs[0]["type"] == "ribo":
                if modif == 0: ribo_binding_reaction(a1_ex, a1, r_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
                elif modif == 1: ribo_binding_reaction(a1_ex, a1, r_1_b, P_in, drugs[0]["P_out"], K_on, K_off, Lambda)
                elif modif == 2:ribo_binding_reaction(a1_ex, a1, r_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on_exp, K_off_exp, Lambda)

        if len(drugs) > 1:
            # こちらは多剤のみ
            K_on_exp = K_on * 1 / (1 + a1 / drugs[0]["K_ma"])
            # K_off_exp = K_on_exp * K_D
            K_off_exp = K_off
            P_in = drugs[1]["P_in"] * 1 / (1 + a1 / drugs[0]["K_ma"])

            if drugs[1]["type"] == "30s":
                if modif == 0: r30_binding_reaction(a2_ex, a2, r30_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on, K_off, Lambda)
                elif modif == 1: r30_binding_reaction(a2_ex, a2, r30_2_b, P_in, drugs[1]["P_out"], K_on, K_off, Lambda)
                elif modif == 2:r30_binding_reaction(a2_ex, a2, r30_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on_exp, K_off_exp, Lambda)
            elif drugs[1]["type"] == "50s":
                if modif == 0: r50_binding_reaction(a2_ex, a2, r50_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on, K_off, Lambda)
                elif modif == 1: r50_binding_reaction(a2_ex, a2, r50_2_b, P_in, drugs[1]["P_out"], K_on, K_off, Lambda)
                elif modif == 2:r50_binding_reaction(a2_ex, a2, r50_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on_exp, K_off_exp, Lambda)
            elif drugs[1]["type"] == "ribo":
                if modif == 0: ribo_binding_reaction(a2_ex, a2, r_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on, K_off, Lambda)
                elif modif == 1: ribo_binding_reaction(a2_ex, a2, r_2_b, P_in, drugs[1]["P_out"], K_on, K_off, Lambda)
                elif modif == 2:ribo_binding_reaction(a2_ex, a2, r_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on_exp, K_off_exp, Lambda)

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

    drugData = {"name": drugName,
                "type": types[drugName],
                "dose": .0,
                "Lambda_0_a": Lambda_0_a[drugName],
                "IC50": IC50[drugName][medium],
                "IC50_a": IC50_a[drugName],
                "K_ma": 15.
                }


    return drugData


def run(drugs=[], step=50., medium=0, legend=[], inpData={}, y0={"r30_u": 30., "r50_u": 30., "r_u": 30., "r_b": .0}):
    """
    ribosomeモデル実行関数

    drugs: [{"name":, "type":, "dose":, }]
    step:
    legend:
    inpData: createModelに渡すDataset
    y0:
    """
    Lambda_0 =  [1.35, 0.85, 0.40]
    dataset = {"r_max": 65.8, "r_min": 19.3, "K_D": 1.0, "K_t": 6.1*10**-2, "K_on": 3.0, "Lambda_0": Lambda_0[medium], "Kd": 1., "p": 1.} # 基本のデータセット
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


# def checkLinerType(inp, eps_rel, pattern, buffpoint=0.):
#     """
#     inp: growth rate list
#     eps_rel: 差の許容率
#     pattern: linertypeを計算するときのパターン
#     """
#     inp = np.array(inp)
#     dinp = inp[1:] - inp[:-1] # 差のリスト
#     before_type = 0 # -1: <, 0: ==, 1: >
#     linertype = 0 # -1: synergistic, 0: additive, 2: antagonistic`
#
#     if pattern == 0:
#         for i in range(len(dinp)):
#             if abs(dinp[i]) > inp[i] * eps_rel: #
#                 if dinp[i] < 0:
#                     current_type = -1 # 減少
#                 else:
#                     current_type = 1 # 増加
#             else:
#                 current_type = 0
#
#             if i == 0:
#                 before_type = current_type
#             else:
#                 if current_type == 0 or before_type == 0 or current_type == before_type:
#                     pass
#                 else:
#                     linerList = np.linspace(inp[0], inp[-1], len(inp)) # 両端点を線形で結んだList
#                     linertype = linerList[i] - inp[i] #
#
#                     return linertype
#                     break
#
#     elif pattern == 1:
#         upper_bound = max(inp[0], inp[-1])
#         lower_bound = min(inp[0], inp[-1])
#         max_inp = max(inp)
#         min_inp = min(inp)
#         if max_inp > upper_bound: # antagonistic
#             linertype = upper_bound - max_inp
#         elif min_inp < lower_bound: # synergistic
#             linertype = lower_bound - min_inp
#
#     elif pattern == 2:
#         if buffpoint != 0:
#             upper_bound = max(inp[0], inp[-1])
#             lower_bound = min(inp[0], inp[-1])
#             max_inp = max(inp)
#             min_inp = min(inp)
#             if max_inp > upper_bound: # antagonistic
#                 linertype = -(max_inp / buffpoint)
#             elif min_inp < lower_bound: # synergistic
#                 linertype = lower_bound - min_inp
#
#     return linertype

def checkLinerType(inp, buff=0.):
    """
    inp: growth rate list
    buff: 差の許容率
    """
    upper_bound = max(inp[0], inp[-1]) + buff
    lower_bound = min(inp[0], inp[-1]) - buff
    max_inp = max(inp)
    min_inp = min(inp)
    linertype = 0
    if max_inp > upper_bound: # antagonistic
        linertype = upper_bound - max_inp
    elif min_inp < lower_bound: # synergistic
        linertype = lower_bound - min_inp
    return linertype


def calcIC(drugNames, a_ex, target, modif=0):
    """
    二分法でIC〜〜を計算
    """
    calc_result = {}
    for dName in drugNames:
        drugs = [makeDrugDatas(dName, 0)] # 薬剤データの作成
        dose_max = a_ex[dName] * 3
        dose_min = .0
        result = 1.
        dose = .0
        while abs(target - abs(result)) > 0.01:
            dose = (dose_max + dose_min) / 2.
            drugs[0]["dose"] = dose
            result, legend = run(drugs, step=100, inpData={"modif": modif}, legend=["r_u"])
            result = calcGrowthrate(result[-1][1])
            if result < target:
                dose_max = dose
            else:
                dose_min = dose
        calc_result[dName] = dose
    return calc_result


def calcBufferingPoint(drugNames, doses):
    """
    二分法でBufferingPointを計算
    """
    a_range = [0, doses[0]]
    a_mid = 0
    b_range = [0, doses[1]]
    b_mid = doses[1]
    result = 1
    drugA = [makeDrugDatas(drugNames[0])]
    drugB = [makeDrugDatas(drugNames[1])]
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


# modules
def makedir(dirname):
    """ディレクトリ作成関数"""
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    del(os)

def makeCmap(c_range={"red": 30, "pink": 5, "white": 10, "light_green": 5, "green": 13, "blue": 17},
            c_num = {"red": "#ff0000", "pink": "#ffc0cb", "white": "#ffffff", "light_green": "#90ee90", "green": "#008000", "blue": "#0000ff"},
            c_list = ["red", "pink", "white", "light_green", "green", "blue"]): # eをかませるためのカラーマップを作る関数
    """自分で定義したカラーマップを返す(固定)"""
    import matplotlib

    c_result = []
    for color in c_list:
        for i in range(c_range[color]):
            c_result.append(c_num[color])
    cmap = matplotlib.colors.ListedColormap(c_result)
    del(matplotlib)
    return cmap

def generate_cmap(colors=["mediumblue", "white", "orangered"]):
    """自分で定義したカラーマップを返す(線形補完)"""
    from matplotlib.colors import LinearSegmentedColormap
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


def checkEpsilon(x, y, xy):
    """
    Nature genetics 2006のepsilonを計算する関数
    e = (Wxy - WxWy)/|min(Wx, Wy) - WxWy|
    """
    result = (xy - x * y) / abs(min(x, y) - x * y)
    return result


def doseResponse(drugs, dose, inpData={"modif": 0}):
    """
        run して，Growth rateを返す関数
    """
    modif2_K_ma = {'Kanamycin': 8.9, 'Streptmycin': 9.2, 'Chloramphenicol': 23.1, 'Tetracycline': 24.7}

    for i in range(len(drugs)):
        drugs[i]["dose"] = dose[i]
        # modifの数値によって，K_maを変更
        if inpData["modif"] < 2: drugs[i]["K_ma"] = 15.
        elif inpData["modif"] == 2: drugs[i]["K_ma"] = modif2_K_ma[drugs[i]["name"]]
    result, legend = run(drugs, step=100, inpData=inpData, legend=["r_u"])
    result = calcGrowthrate(result[-1][1])
    return result

def createSlopedose(slope, midPointList, divnum=11):
    """
        midPointに応じて，該当するSlopeで作成したDoseの組み合わせリストを返す関数
    """
    doseX = np.linspace(0, midPointList[0] * (1 + slope), divnum)
    doseY = np.linspace(0, midPointList[1] * (1 + (1 / slope)), divnum)[::-1]
    return_list = [[doseX[i], doseY[i]] for i in range(len(doseX))]
    return [[doseX[i], doseY[i]] for i in range(len(doseX))]

def calcEpsilon(drugNames, doses, inpData={"modif": 1}):
    """
    """
    result_list = []
    for i in range(3):
        if i < 2:
            drugs = [makeDrugDatas(drugNames[i])]
            result_list.append(doseResponse(drugs, [doses[i]], inpData=inpData))
        else:
            drugs = [makeDrugDatas(drugNames[0]), makeDrugDatas(drugNames[1])]
            result_list.append(doseResponse(drugs, doses, inpData=inpData))

    return epsilon(result_list[0], result_list[1], result_list[2])

def calcGrowthRate(a, r_min=19.3, K_t=6.1*10**-2, Lambda_0=1.35):
    """
    growth_rateを計算する関数
    """
    result = (a - r_min) * K_t / Lambda_0
    return result

def sim(drugs, dose, inpData={}):
    for i in range(len(drugs)):
        drugs[i]["dose"] = dose[i]
    result, legend = run(drugs, legend=["r_u"], inpData=inpData)
    return calcGrowthRate(result[-1][1])

def setTickLabel(data, ax):
    a1DoseList = list(set(data["a1"].tolist()))[::2] # y軸に表示したい数のリスト
    a2DoseList = list(set(data["a2"].tolist()))[::2] # x軸に表示したい数のリスト

    ax.set_xticks(list(np.linspace(0.5, 10.5, len(a2DoseList)))) # xticksのせてい
    ax.set_xticklabels(list(map(str, a2DoseList)))
    ax.set_yticks(list(np.linspace(0.5, 10.5, len(a1DoseList))))
    ax.set_yticklabels(list(map(str, a1DoseList)))

def divideDoses(drugNames, IC30, num, length=101, splitNum=101):
    """
        doseのリストを分割加工する関数
        drugNames : 薬剤の名前のリスト
        IC30      : 
        num       : 何番目のリストか
        length    : 何間隔でやっているか
        splitNum  : 全部で何 * 何のシミュレーションになるか
    """
    doses = [[x, y] 
             for x in np.linspace(0, IC30[drugNames[0]] * 2, splitNum) 
             for y in np.linspace(0, IC30[drugNames[1]] * 2, splitNum)]
    if (num + 1) * length < len(doses):
        doses = doses[num * length : (num + 1) * length]
    else :
        doses = doses[num * length :]

    return doses

def divideDosesNeweval(drugNames, IC30, num, length=5, splitNum=11):
    """
        newEvalのdoseリストを分割加工する関数
        drugNames : 薬剤の名前のリスト
        IC30      :
        num       : 何番目のリストか
        length    : 何間隔でやっているか
        splitNum  : 何分割でLinerTypeを作成するか
    """
    slopeList = [1./4, 1./2., 1., 2., 4.]
    midPointList = [x for x in itr.zip_longest(np.linspace(0, IC30[drugName[0]], splitNum), np.linspace(0, IC30[drugName[1]], splitNum))][1:] # 0, 0を除いた10点．中点なので，IC30でOK．
    doseList = [[[x for x in itr.zip_longest(np.linspace(0, midPoint[0] * (1 + slope), splitNum), np.linspace(0, midPoint[1] * (1 + (1 / slope)), splitNum)[::-1])] for midPoint in midPointList] for slope in slopeList]
    doses = doseList[num]

    return doses

def sim_comb(drugs, doses, inpData={"modif": 0}, target=None):
    if target:
        drugs[0]["type"] = target[0]
        drugs[1]["type"] = target[1]
    resultList = []
    for index, dose in enumerate(doses):
        print("    step: {} >> ".format(index))
        resultList.append([dose[0], dose[1], sim(drugs, dose, inpData)])
    data = pd.DataFrame(resultList, columns=["a1", "a2", "growth"])
    return data

def sim_oldeval(drugName, doses, inpData={"modif": 0}, target=None):
    """
        論文の評価関数を使ったシミュレーション．
        drugName : 薬剤の名前のリスト
        doses    :
        target   : 
    """
    drug1 = [makeDrugDatas(drugName[0])]
    drug2 = [makeDrugDatas(drugName[1])]
    drug3 = [makeDrugDatas(drugName[0]), makeDrugDatas(drugName[1])]
    if target:
        drug1[0]["type"] = target[0]
        drug2[0]["type"] = target[1]
        drug3[0]["type"] = target[0]
        drug3[1]["type"] = target[1]
    resultList = []
    for index, dose in enumerate(doses):
        print("    step: {} >> ".format(index))
        x = sim(drug1, [dose[0]], inpData)
        y = sim(drug2, [dose[1]], inpData)
        xy = sim(drug3, dose, inpData)
        if 0 in dose: Eps = 0
        else: Eps = checkEpsilon(x, y, xy)
        resultList.append([dose[0], dose[1], x, y, xy, Eps])
    data = pd.DataFrame(resultList, columns=["a1", "a2", "growth1", "growth2", "growth3", "epsilon"])
    return data

def sim_neweval(drug, slope, doseList, inpData={"modif": 0}, target=None):
    if target:
        for i in range(len(target)):
            drug[i]["type"] = target[i]
    resultList = []
    for index, doses in enumerate(doseList):
        print("    step : {} >> ".format(index)) 
        growthList = [sim(drugs, dose, inpData) for dose in doses]
        resultList.append([(index + 1) * 10, slope, growthList, checkLinerType(growthList)])
    data = pd.DataFrame(resultList, columns=["MidPoint", "Slope", "resultList", "LinerType"])
    return data


if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]

    # IC30の計算
    makedir("IC30")
    IC30_file_modif1 = "IC30/ribo5_modif1.csv"
    IC30 = [] # 0: modif1, 1: modif2
    try:
        IC30_df = pd.read_csv(IC30_file_modif1)
        IC30.append({i: IC30_df[i][0] for i in IC30_df.columns})
    except:
        IC30.append(calcIC(drugNames, {drugName: 20 for drugName in drugNames}, .3, modif=1))
        IC30_df = pd.DataFrame({i: [IC30[0][i]] for i in drugNames})
        IC30_df.to_csv(IC30_file_modif1, index=False)
    
    IC30_file_modif2 = "IC30/ribo5_modif2.csv"
    try:
        IC30_df = pd.read_csv(IC30_file_modif2)
        IC30.append({i: IC30_df[i][0] for i in IC30_df.columns})
    except:
        IC30.append(calcIC(drugNames, {drugName: 20 for drugName in drugNames}, .3, modif=2))
        IC30_df = pd.DataFrame({i: [IC30[1][i]] for i in drugNames})
        IC30_df.to_csv(IC30_file_modif2, index=False)
    
    # ヒートマップ用データ作成
    # csvdir = "results/ribo5/csv/sim100"
    # makedir(csvdir)
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # num = int(sys.argv[-1])
    # print("start simulation >> ")
    # for modif in range(1, 3):
    #     print("  modif {} >> ".format(modif))
    #     for drugName in drugNameList:
    #         print("    {} vs {} >> ".format(drugName[0], drugName[1]))
    #         inpData = {"modif": modif}
    #         dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(drugName))
    #         makedir(dirName)
    #         drugs = [makeDrugDatas(drugName[0]), makeDrugDatas(drugName[1])]
    #         doses = divideDoses(drugName, IC30[modif - 1], num, 101, 101)
    #         df = sim_comb(drugs, doses, inpData)
    #         df.to_csv("{}/{}.csv".format(dirName, num), index=False)

    # ヒートマップ用データ作成（仮想薬剤）
    # csvdir = "results/ribo5/csv/sim100_v"
    # makedir(csvdir)
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["30s", "30s"], ["30s", "50s"]]
    # num = int(sys.argv[-1])
    # print("start combination >> ")
    # for modif in range(1, 3):
    #     print("  modif {} >>".format(modif))
    #     inpData = {"modif": modif}
    #     for drugName in drugNameList:
    #         print("    {} vs {} >> ".format(drugName[0], drugName[1]))
    #         doses = divideDoses(drugName, IC30[modif - 1], num, 101, 101)
    #         for target in targetList:
    #             dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
    #             makedir(dirName)
    #             print("      {} vs {} >> ".format(target[0], target[1]))
    #             drugs = [makeDrugDatas(drugName[0]), makeDrugDatas(drugName[1])]
    #             df = sim_comb(drugs, doses, inpData, target)
    #             df.to_csv("{}/{}.csv".format(dirName, num), index = False)
    
    # simulate oldeval
    # csvdir = "results/ribo5/csv/old100"
    # makedir(csvdir)
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # num = int(sys.argv[-1])
    # print("start simulation >> ")
    # for modif in range(1, 3):
    #     print("  modif {} >> ".format(modif))
    #     for drugName in drugNameList:
    #         print("    {} vs {} >> ".format(drugName[0], drugName[1]))
    #         inpData = {"modif": modif}
    #         dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(drugName))
    #         makedir(dirName)
    #         doses = divideDoses(drugName, IC30[modif - 1], num, 101, 101)
    #         df = sim_oldeval(drugName, doses, inpData)
    #         df.to_csv("{}/{}.csv".format(dirName, num), index=False)

    # simulate oldeval (virtual drug)
    # csvdir = "results/ribo5/csv/old100_v"
    # makedir(csvdir)
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["30s", "30s"], ["30s", "50s"]]
    # num = int(sys.argv[-1])
    # print("start combination >> ")
    # for modif in range(1, 3):
    #     print("  modif {} >>".format(modif))
    #     inpData = {"modif": modif}
    #     for drugName in drugNameList:
    #         print("    {} vs {} >> ".format(drugName[0], drugName[1]))
    #         doses = divideDoses(drugName, IC30[modif - 1], num, 101, 101)
    #         for target in targetList:
    #             dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
    #             makedir(dirName)
    #             print("      {} vs {} >> ".format(target[0], target[1]))
    #             df = sim_oldeval(drugName, doses, inpData, target)
    #             df.to_csv("{}/{}.csv".format(dirName, num), index = False)

    ## neweval simulation
    csvdir = "results/ribo5/csv/new100"
    makedir(csvdir)
    num = int(sys.argv[-1])
    drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    slopeList = [1./4., 1./2., 1., 2., 4.]
    print("start simulation >> ")
    for modif in range(1, 3):
        print("  modif {} >> ".format(modif))
        inpData = {"modif": modif}
        for drugName in drugNameList:
            print("    {} vs {} >>".format(drugName[0], drugName[1]))
            dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(drugName))
            makedir(dirName)
            doses = divideDosesNeweval(drugName, IC30[modif - 1], num, 5, 11)
            drugs = [makeDrugDatas(drugName[0]), makeDrugDatas(drugName[1])]
            df = sim_neweval(drugs, slopeList[num], doses, inpData)
            df.to_csv("{}/{}.csv".format(dirName, num), index=False)
    
    ## neweval simulation (virtual drug)
    csvdir = "results/ribo5/csv/new100_v"
    makedir(csvdir)
    num = int(sys.argv[-1])
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["30s", "30s"], ["30s", "50s"]]
    slopeList = [1./4., 1./2., 1., 2., 4.]
    print("start combination >> ")
    for modif in range(1, 3):
        print("  modif {} >> ".format(modif))
        inpData = {"modif": modif}
        for drugName in drugNameList:
            print("  {} vs {} >> ".format(drugName[0], drugName[1]))
            doses = divideDosesNeweval(drugName, IC30[modif - 1], num)
            drugs = [makeDrugDatas(drugName[0]), makeDrugDatas(drugName[1])]
            for target in targetList:
                print("    {} vs {} >> ".format(target[0], target[1]))
                dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
                makedir(dirName)
                df = sim_neweval(drugs, slopeList[num], doses, inpData, target)
                df.to_csv("{}/{}.csv".format(dirName, num), index=False)
    

