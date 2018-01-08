# coding: utf-8
"""
    r30, r50に，それぞれ2つ薬剤が結合できるようにしたもの．
    State : no_reaction
"""

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True
import itertools as itr
import sys
from ribo4 import makedir
from ribo4 import checkEpsilon
from ribo4 import checkLinerType

# main
@reaction_rules
def r30_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda, r_min, num, addReaction):
    ~a_ex == a | (P_in * a_ex, P_out)
    if num == 0: # 一つ目の薬剤の場合
        a + r30_u_u == r30_b_u | (K_on, K_off)
        a + r30_u_b == r30_b_b | (K_on, K_off)
        if addReaction: a + r_u > r30_b_u + r50_u_u | K_on * a * (r_u - r_min) # add
        r30_b_u > ~r30_b_u | Lambda * r30_b_u
    elif num == 1: # 二つ目の薬剤の場合
        a + r30_u_u == r30_u_b | (K_on, K_off)
        a + r30_b_u == r30_b_b | (K_on, K_off)
        if addReaction: a + r_u > r30_u_b + r50_u_u | K_on * a * (r_u - r_min) # add
        r30_u_b > ~r30_u_b | Lambda * r30_u_b
    a > ~a | Lambda * a

@reaction_rules
def r50_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda, r_min, num, addReaction):
    ~a_ex == a | (P_in * a_ex, P_out)
    if num == 0: # 一つ目の薬剤の場合
        a + r50_u_u == r50_b_u | (K_on, K_off)
        a + r50_u_b == r50_b_b | (K_on, K_off)
        if addReaction: a + r_u > r50_b_u + r30_u_u | K_on * a * (r_u - r_min) # add
        r50_b_u > ~r50_b_u | Lambda * r50_b_u
    elif num == 1: # 二つ目の薬剤の場合
        a + r50_u_u == r50_u_b | (K_on, K_off)
        a + r50_b_u == r50_b_b | (K_on, K_off)
        if addReaction: a + r_u > r50_u_b + r30_u_u | K_on * a * (r_u - r_min) # add
        r50_u_b > ~r50_u_b | Lambda * r50_u_b
    a > ~a | Lambda * a

@reaction_rules
def ribo_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda, r_min, num, addReaction):
    ~a_ex == a | (P_in * a_ex, P_out)
    a + r_u == r_b | (K_on * a * (r_u - r_min), K_off)
    r_b > a + r30_u_u + r50_u_u | Kd
    a > ~a | Lambda * a

def createModel(drugs=[], dataset={}, addReaction=False):
    # 定数
    K_t = 6.1 * 10 ** -2

    # 薬剤関連
    K_D = 1. # 薬剤の解離定数と結合定数の比
    K_on = 3.0 # 薬剤の結合定数
    K_off = K_on * K_D # 薬剤の解離定数

    # リボソーム，growth関連
    medium = 0 # 培地の設定 0: Gly, 1: Gly_CAA, 2: Gly_RDM
    Lambda_0_List =  [1.35, 0.85, 0.40] # 初期培地でのGrowth Rateのリスト．内容はmediumに準拠．
    # mediumに応じて初期培地でのGrowth Rateを変更
    Lambda_0 = Lambda_0_List[medium] # 初期培地でのGrowth Rate
    r_min = 19.3 # リボソームの最小値
    r_max = 65.8 # リボソームの最大値
    Delta_r = r_max - r_min # リボソームの最大値と最小値の差
    r_u_0 = Lambda_0 / K_t + r_min # 定常状態でのr_uの値の予測

    # リボソームサブユニット関連
    Kd = 1. # リボソームサブユニットの解離定数
    p = 1. # リボソームサブユニットの存在する比率(r_uと比べて)
    Ka = (Kd / K_t + r_u_0) * Lambda_0 / ((p * r_u_0) ** 2) # リボソームサブユニットの結合定数

    # 引数でdatasetを受け取った場合，更新する
    if dataset:
        for key, value in dataset.items():
            exec("{} = {}".format(key, value))

    with reaction_rules():
        # expression
        Lambda = (r_u - r_min) * K_t # Growth rate
        SUP = (Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1 / K_t / Delta_r))) * (1 + p) # subunit product expression

        # drug reaction
        if drugs:
            for index, drug in enumerate(drugs):
                a_ex = _eval("a{}_ex".format(index))
                a = _eval("a{}".format(index))
                P_in = (r_max - r_min) * drug["Lambda_0_a"] / 2.0 / drug["IC50_a"] # 薬剤の流入
                P_out = (drug["Lambda_0_a"] / 2) ** 2.0 / K_t / K_D # 薬剤の流出
                exec('{}_binding_reaction(a_ex, a, P_in, P_out,\
                                           K_on, K_off, Lambda, r_min, index, addReaction)'\
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
    """
        drug: drug data (Data after having changed makeDrugDatas)
    """
    model = createModel(drugs, inpData)

    # check reaction rules
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

def sim(drugs, dose):
    # runなどのシミュレーションをして，結果を返す関数
    for i in range(len(drugs)):
        drugs[i]["dose"] = dose[i]
    result, legend = run(drugs, legend=["r_u"])
    return calcGrowthRate(result[-1][1])

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
                "K_ma": 15. # antagonisticにした際のフィッティングパラメータ
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


def sim_comb(drugs, doses, target=None):
    if target:
        drugs[0]["target"] = target[0]
        drugs[1]["target"] = target[1]
    resultList = []
    for index, dose in enumerate(doses):
        print("    step: {} >> ".format(index))
        resultList.append([dose[0], dose[1], sim(drugs, dose)])
    data = pd.DataFrame(resultList, columns=["a1", "a2", "growth"])
    return data

def sim_oldeval(drugName, doses, target=None):
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
        drug1[0]["target"] = target[0]
        drug2[0]["target"] = target[1]
        drug3[0]["target"] = target[0]
        drug3[1]["target"] = target[1]
    resultList = []
    for index, dose in enumerate(doses):
        print("    step: {} >> ".format(index))
        x = sim(drug1, [dose[0]])
        y = sim(drug2, [dose[1]])
        xy = sim(drug3, dose)
        if 0 in dose: Eps = 0
        else: Eps = checkEpsilon(x, y, xy)
        resultList.append([dose[0], dose[1], x, y, xy, Eps])
    data = pd.DataFrame(resultList, columns=["a1", "a2", "growth1", "growth2", "growth3", "epsilon"])
    return data

def sim_neweval(drug, slope, doseList, target=None):
    if target:
        for i in range(len(target)):
            drug[i]["target"] = target[i]
    resultList = []
    for index, doses in enumerate(doseList):
        growthList = [sim(drugs, dose) for dose in doses]
        resultList.append([(index + 1) * 10, slope, growthList, checkLinerType(growthList)])
    data = pd.DataFrame(resultList, columns=["MidPoint", "Slope", "resultList", "LinerType"])
    return data

if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]

    # IC30の計算
    makedir("IC30")
    IC30_file = "IC30/ribo7.csv"
    try:
        IC30_df = pd.read_csv(IC30_file)
        IC30 = {i: IC30_df[i][0] for i in IC30_df.columns}
    except:
        IC30 = calcIC(drugNames, {drugName: 100 for drugName in drugNames}, .3)
        IC30_df = pd.DataFrame({i: [IC30[i]] for i in drugNames})
        IC30_df.to_csv(IC30_file, index=False)
    
    ## double simulation
    # csvdir = "./results/ribo7/csv/sim100"
    # makedir(csvdir)
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # num = int(sys.argv[-1])
    # print("start combination >> ")
    # for drugName in drugNameList:
    #     dirName = "{}/{}".format(csvdir, "_".join(drugName))
    #     makedir(dirName)
    #     print("{} vs {}".format(drugName[0], drugName[1]))
    #     drugs = [makeDrugDatas(drugName[0]), makeDrugDatas(drugName[1])]
    #     doses = divideDoses(drugName, IC30, num, 101, 101)
    #     df = sim_comb(drugs, doses)
    #     df.to_csv("{}/{}_{}.csv".format(dirName, "_".join(drugName), num), index=False)
    
    ## double simulation (virtual drugs)
    # csvdir = "./results/ribo7/csv/sim100_v"
    # makedir(csvdir)
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["r30", "r30"], ["r30", "r50"]]
    # num = int(sys.argv[-1])
    # print("start combination >> ")
    # for drugName in drugNameList:
    #     print("  {} vs {} >> ".format(drugName[0], drugName[1]))
    #     doses = divideDoses(drugName, IC30, num, 101, 101)
    #     for target in targetList:
    #         dirName = "{}/{}".format(csvdir, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
    #         makedir(dirName)
    #         print("    {} vs {} >> ".format(target[0], target[1]))
    #         drugs = [makeDrugDatas(drugName[0]), makeDrugDatas(drugName[1])]
    #         df = sim_comb(drugs, doses, target)
    #         df.to_csv("{}/{}.csv".format(dirName, num), index = False)

    ## oldeval simulation
    # csvdir = "./results/ribo7/csv/old100"
    # makedir(csvdir)
    # num = int(sys.argv[-1])
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # print("start simulation >> ")
    # for drugName in drugNameList:
    #     print("  {} vs {}".format(drugName[0], drugName[1]))
    #     dirName = "{}/{}".format(csvdir, "_".join(drugName))
    #     makedir(dirName)
    #     doses = divideDoses(drugName, IC30, num, 101, 101)
    #     df = sim_oldeval(drugName, doses)
    #     df.to_csv("{}/{}.csv".format(dirName, num), index=False)
        
    ## oldeval simulation (virtual drug)
    # csvdir = "results/ribo7/csv/old100_v"
    # makedir(csvdir)
    # num = int(sys.argv[-1])
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["r30", "r30"], ["r30", "r50"]]
    # print("start combination >> ")
    # for drugName in drugNameList:
    #     print("  {} vs {} >> ".format(drugName[0], drugName[1]))
    #     doses = divideDoses(drugName, IC30, num, 101, 101)
    #     for target in targetList:
    #         dirName = "{}/{}".format(csvdir, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
    #         makedir(dirName)
    #         print("    {} vs {} >> ".format(target[0], target[1]))
    #         df = sim_oldeval(drugName, doses, target)
    #         df.to_csv("{}/{}.csv".format(dirName, num), index = False)
    
    ## neweval simulation
    # csvdir = "results/ribo7/csv/new100"
    # makedir(csvdir)
    # num = int(sys.argv[-1])
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # slopeList = [1./4., 1./2., 1., 2., 4.]
    # print("start simulation >> ")
    # for drugName in drugNameList:
    #     dirName = "{}/{}".format(csvdir, "_".join(drugName))
    #     makedir(dirName)
    #     print("  {} vs {} >>".format(drugName[0], drugName[1]))
    #     doses = divideDosesNeweval(drugName, IC30, num, 5, 11)
    #     drugs = [makeDrugDatas(drugName[0]), makeDrugDatas(drugName[1])]
    #     df = sim_neweval(drugs, slopeList[num], doses)
    #     df.to_csv("{}/{}.csv".format(dirName, num), index=False)

    ## neweval simulation (virtual drug)
    csvdir = "results/ribo7/csv/new100_v"
    makedir(csvdir)
    num = int(sys.argv[-1])
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["r30", "r30"], ["r30", "r50"]]
    slopeList = [1./4., 1./2., 1., 2., 4.]
    print("start combination >> ")
    for drugName in drugNameList:
        print("  {} vs {} >> ".format(drugName[0], drugName[1]))
        doses = divideDosesNeweval(drugName, IC30, num)
        drugs = [makeDrugDatas(drugName[0]), makeDrugDatas(drugName[1])]
        for target in targetList:
            print("    {} vs {} >> ".format(target[0], target[1]))
            dirName = "{}/{}".format(csvdir, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
            makedir(dirName)
            df = sim_neweval(drugs, slopeList[num], doses, target)
            df.to_csv("{}/{}.csv".format(dirName, num), index=False)

