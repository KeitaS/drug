# coding: utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
import math
import itertools as itr
import seaborn as sns
import pandas as pd
import sys
util.decorator.SEAMLESS_RATELAW_SUPPORT = True

@reaction_rules
def createVariable(drug_ex, drug, P_in, P_out):
    # 薬剤の流入式を入れる関数
    ~drug_ex == drug | (P_in * drug_ex, P_out * drug)

def createDrugData(drugName):
    # 薬剤のプロファイルを作成する関数
    # すべて，Lambda_0 = 1.35
    datas = {
        "Lambda_0_a" : {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83},
        "IC50_a"     : {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49},
        "target"     : {"Streptmycin": "A", "Kanamycin": "A", "Tetracycline": "A", "Chloramphenicol": "C"}
    }

    drugData = {
        "name"       : drugName,
        "target"     : datas["target"][drugName],
        "dose"       : .0,
        "Lambda_0_a" : datas["Lambda_0_a"][drugName],
        "IC50_a"     : datas["IC50_a"][drugName]
    }
    return drugData

def createModel(drugs, K_D=1., K_on=3., sub_k_d=1., sub_p=1., Lambda_0=1.35):
    """
    モデルを作成する関数
    coding rules
    薬剤に関連する変数はアッパースネーク，それ以外は基本的にローワースネーク．
    一部命名規則の問題から，アッパースネークを利用．
    r_max: リボソームの最大収容量(µM)．
    r_min: リボソームの最小収容量(µM)．
    r_u_0: リボソームの初期値．
    sub_k_a: サブユニットの結合定数．
    sub_k_d: サブユニットの解離定数．
    sub_p: サブユニットのリボソームに対する存在比率．
    K_D: 薬剤の平衡解離定数．
    K_t: 翻訳能力．empirical relations，
    K_on: 薬剤の結合定数．
    K_off: 薬剤の解離定数．
    P_in: 薬剤の流入速度．
    P_off: 薬剤の排出速度．
    """
    # 定数
    K_t = 6.1 * 10 ** -2 # 固定
    r_min = 19.3 # 固定 元論文から定義の変更
    r_max = 65.8 # 固定
    Lambda_0_list = [1.35, 0.85, 0.4]

    delta_r = r_max - (r_min * (1 + sub_p)) # 元論文から定義の変更
    K_off = K_on * K_D
    r_u_0 = Lambda_0 / K_t + r_min
    sub_k_a = (sub_k_d / K_t + r_u_0) * Lambda_0 / ((sub_p * r_u_0) ** 2)


    with reaction_rules():
        r30_tot = (r30   + r30_r50   + r30_r50C   + r30_r50D   + r30_r50CD +
                 r30A  + r30A_r50  + r30A_r50C  + r30A_r50D  + r30A_r50CD +
                 r30B  + r30B_r50  + r30B_r50C  + r30B_r50D  + r30B_r50CD +
                 r30AB + r30AB_r50 + r30AB_r50C + r30AB_r50D + r30AB_r50CD)
        Lambda = (r30_tot - r_min) * K_t * r30_r50 / r30_tot
        SUP = Lambda * r30_tot
        # SUP = delta_r * ((1 + sub_p) / (K_t * delta_r) - 1 / Lambda_0)
        for index, drug in enumerate(drugs):
            # 薬剤の流入の式を追加
            drug_ex = _eval("drug{}_ex".format(index))
            target = _eval(drug["target"])
            P_in = delta_r * drug["Lambda_0_a"] / 2.0 / drug["IC50_a"] # 薬剤の流入
            P_out = (drug["Lambda_0_a"] / 2) ** 2.0 / K_t / K_D # 薬剤の流出
            createVariable(drug_ex, target, P_in, P_out)

        # DrugA(30s)
        # 薬剤希釈
        A > ~A | Lambda * A
        # 薬剤結合
        ## 標的サブユニットへの結合
        A + r30  == r30A  | (K_on, K_off)
        A + r30B == r30AB | (K_on, K_off)
        ## 機能性リボソームへの結合
        A + r30_r50    == r30A_r50    | (K_on, K_off)
        A + r30_r50C   == r30A_r50C   | (K_on, K_off)
        A + r30_r50D   == r30A_r50D   | (K_on, K_off)
        A + r30_r50CD  == r30A_r50CD  | (K_on, K_off)
        A + r30B_r50   == r30AB_r50   | (K_on, K_off)
        A + r30B_r50C  == r30AB_r50C  | (K_on, K_off)
        A + r30B_r50D  == r30AB_r50D  | (K_on, K_off)
        A + r30B_r50CD == r30AB_r50CD | (K_on, K_off)

        # DrugB(30s)
        # 薬剤希釈
        B > ~B | Lambda * B
        # 薬剤結合
        ## 標的サブユニットへの結合
        B + r30  == r30B  | (K_on, K_off)
        B + r30A == r30AB | (K_on, K_off)
        ## 機能性リボソームへの結合
        B + r30_r50    == r30B_r50    | (K_on, K_off)
        B + r30_r50C   == r30B_r50C   | (K_on, K_off)
        B + r30_r50D   == r30B_r50D   | (K_on, K_off)
        B + r30_r50CD  == r30B_r50CD  | (K_on, K_off)
        B + r30A_r50   == r30AB_r50   | (K_on, K_off)
        B + r30A_r50C  == r30AB_r50C  | (K_on, K_off)
        B + r30A_r50D  == r30AB_r50D  | (K_on, K_off)
        B + r30A_r50CD == r30AB_r50CD | (K_on, K_off)

        # DrugC(50s)
        # 薬剤希釈
        C > ~C | Lambda * C
        # 薬剤結合
        ## 標的サブユニットへの結合
        C + r50  == r50C  | (K_on, K_off)
        C + r50D == r50CD | (K_on, K_off)
        ## 機能性リボソームへの結合
        C + r30_r50    == r30_r50C    | (K_on, K_off)
        C + r30_r50D   == r30_r50CD   | (K_on, K_off)
        C + r30A_r50   == r30A_r50C   | (K_on, K_off)
        C + r30A_r50D  == r30A_r50CD  | (K_on, K_off)
        C + r30B_r50   == r30B_r50C   | (K_on, K_off)
        C + r30B_r50D  == r30B_r50CD  | (K_on, K_off)
        C + r30AB_r50  == r30AB_r50C  | (K_on, K_off)
        C + r30AB_r50D == r30AB_r50CD | (K_on, K_off)

        # DrugD(50s)
        # 薬剤希釈
        D > ~D | Lambda * D
        # 薬剤結合
        ## 標的サブユニットへの結合
        D + r50  == r50D  | (K_on, K_off)
        D + r50C == r50CD | (K_on, K_off)
        ## 機能性リボソームへの結合
        D + r30_r50    == r30_r50D    | (K_on, K_off)
        D + r30_r50C   == r30_r50CD   | (K_on, K_off)
        D + r30A_r50   == r30A_r50D   | (K_on, K_off)
        D + r30A_r50C  == r30A_r50CD  | (K_on, K_off)
        D + r30B_r50   == r30B_r50D   | (K_on, K_off)
        D + r30B_r50C  == r30B_r50CD  | (K_on, K_off)
        D + r30AB_r50  == r30AB_r50D  | (K_on, K_off)
        D + r30AB_r50C == r30AB_r50CD | (K_on, K_off)

        # リボソームの結合パターン
        r30   + r50   == r30_r50     | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30_r50     / r30_tot)) # もともとはこれだけ
        r30   + r50C  == r30_r50C    | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30_r50C    / r30_tot))
        r30   + r50D  == r30_r50D    | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30_r50D    / r30_tot))
        r30   + r50CD == r30_r50CD   | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30_r50CD   / r30_tot))
        r30A  + r50   == r30A_r50    | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30A_r50    / r30_tot))
        r30A  + r50C  == r30A_r50C   | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30A_r50C   / r30_tot))
        r30A  + r50D  == r30A_r50D   | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30A_r50D   / r30_tot))
        r30A  + r50CD == r30A_r50CD  | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30A_r50CD  / r30_tot))
        r30B  + r50   == r30B_r50    | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30B_r50    / r30_tot))
        r30B  + r50C  == r30B_r50C   | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30B_r50C   / r30_tot))
        r30B  + r50D  == r30B_r50D   | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30B_r50D   / r30_tot))
        r30B  + r50CD == r30B_r50CD  | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30B_r50CD  / r30_tot))
        r30AB + r50   == r30AB_r50   | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30AB_r50   / r30_tot))
        r30AB + r50C  == r30AB_r50C  | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30AB_r50C  / r30_tot))
        r30AB + r50D  == r30AB_r50D  | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30AB_r50D  / r30_tot))
        r30AB + r50CD == r30AB_r50CD | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30AB_r50CD / r30_tot))

        # リボソームサブユニットの希釈
        r30         > ~r30         | Lambda * r30
        r50         > ~r50         | Lambda * r50
        r30A        > ~r30A        | Lambda * r30A
        r30B        > ~r30B        | Lambda * r30B
        r30AB       > ~r30AB       | Lambda * r30AB
        r50C        > ~r50C        | Lambda * r50C
        r50D        > ~r50D        | Lambda * r50D
        r50CD       > ~r50CD       | Lambda * r50CD
        r30_r50     > ~r30_r50     | Lambda * r30_r50
        r30_r50C    > ~r30_r50C    | Lambda * r30_r50C
        r30_r50D    > ~r30_r50D    | Lambda * r30_r50D
        r30_r50CD   > ~r30_r50CD   | Lambda * r30_r50CD
        r30A_r50    > ~r30A_r50    | Lambda * r30A_r50
        r30A_r50C   > ~r30A_r50C   | Lambda * r30A_r50C
        r30A_r50D   > ~r30A_r50D   | Lambda * r30A_r50D
        r30A_r50CD  > ~r30A_r50CD  | Lambda * r30A_r50CD
        r30B_r50    > ~r30B_r50    | Lambda * r30B_r50
        r30B_r50C   > ~r30B_r50C   | Lambda * r30B_r50C
        r30B_r50D   > ~r30B_r50D   | Lambda * r30B_r50D
        r30B_r50CD  > ~r30B_r50CD  | Lambda * r30B_r50CD
        r30AB_r50   > ~r30AB_r50   | Lambda * r30AB_r50
        r30AB_r50C  > ~r30AB_r50C  | Lambda * r30AB_r50C
        r30AB_r50D  > ~r30AB_r50D  | Lambda * r30AB_r50D
        r30AB_r50CD > ~r30AB_r50CD | Lambda * r30AB_r50CD

        # リボソームサブユニットの生成
        ~r30 > r30 | SUP
        ~r50 > r50 | SUP

    return get_model()

def run(drugs=[], step=50., inpData={}, y0={"r30": 30., "r50": 30., "r30_r50": 30.}, sp_list=["r30_r50"]):
    """
        シミュレーションを実行する関数．
        drugs: 薬剤データがはいったリストを入れる．リストの長さは2まで対応．doseも入れておく．
        step: シミュレーションステップ数．
        inpData: 定数を変えるときに使用する．詳しくは，createModelを参照．
        y0: 特定の値の初期値．デフォルトの値に特に意味はない．
        sp_list: returnでほしいSpeciesをリスト化して渡すために必要．
    """
    inpData["drugs"] = drugs
    model = createModel(**inpData)

    for index, drug in enumerate(drugs):
        y0["drug{}_ex".format(index)] = drug["dose"]

    # for i, rr in enumerate(model.reaction_rules()):
    #     print("{:03d}: {}".format(i, rr.as_string()))
    runsim = run_simulation(
        step,
        solver = "ode",
        y0 = y0,
        return_type = "observer",
        model = model,
        species_list = sp_list
    )

    data = runsim.data()
    return data

def calcGrowthRate(r30_r50, r30_tot, Lambda_0=1.35):
    """
        r30_totとr30_r50を用いてGrowth Rateを計算する関数．
        r_min，K_tは定数なので，入力しない．
    """
    r_min = 19.3
    K_t = 6.1 * 10 ** -2
    Lambda = (r30_tot - r_min) * K_t * (r30_r50 / r30_tot)
    growth = Lambda / Lambda_0
    return growth

def calcGrowth(r30_r50, r30_tot):
    """
        r30_totとr30_r50を用いてGrowth Rateを計算する関数．
        r_min，K_tは定数なので，入力しない．
    """
    r_min = 19.3
    K_t = 6.1 * 10 ** -2
    Lambda = (r30_tot - r_min) * K_t * (r30_r50 / r30_tot)
    return Lambda

def checkEpsilon(x, y, xy): # nature genesis 2006's evaluation
    return (xy - x * y) / abs(min(x, y) - x * y)

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

def sim(drugs, dose, step=50., inpData={}, y0={"r30": 30., "r50": 30., "r30_r50": 30.}, Lambda_0=1.35):
    """
        Runして，Growth Rateを返す関数.
        sp_listはr30_totを計算できるようにリストを作成している．
        データに不備がないよう，すべてのデータをDataFrame型で返すようにしている．
    """
    # doseを入れる
    for i in range(len(drugs)):
        drugs[i]["dose"] = dose[i]
    
    # rtot を計算するために sp_listを作成
    sp_list = ["r30", "r30_r50", "r30_r50C", "r30_r50D", "r30_r50CD",
        "r30A", "r30A_r50", "r30A_r50C", "r30A_r50D", "r30A_r50CD",
        "r30B", "r30B_r50", "r30B_r50C", "r30B_r50D", "r30B_r50CD",
        "r30AB", "r30AB_r50", "r30AB_r50C", "r30AB_r50D", "r30AB_r50CD"]
    # runを使って，tを除去したリストを作成
    result = run(drugs, step, inpData, y0, sp_list=sp_list)[-1][1:]
    r30_r50 = result[1]
    r30_tot = sum(result)
    # growthを計算(growth rateは後々計算)
    growth = calcGrowth(r30_r50, r30_tot)
    # resultにr30_totとgrowthを追加
    result += [r30_tot, growth]
    # columnsにするsp_listにr30_totとgrwthを追加
    sp_list += ["r30_tot", "growth"]
    # DataFrame型でt以外のすべてのデータを返す
    df = pd.DataFrame([result], columns=sp_list)
    return (growth, df)

def calcIC(dNames, a_ex, target):
    """
    二分法でIC〜〜を計算
    """
    calc_result = {}
    for dName in dNames:
        print(dName)
        drugs = [createDrugData(dName)] # 薬剤データの作成
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
            result = sim(drugs)[0]
            if result < target:
                dose_max = dose
            else:
                dose_min = dose
        calc_result[dName] = dose
    return calc_result

def makedir(dirname):
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    del(os)

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

def sim_comb(drugs, doses, Lambda_0, target=None):
    """
        多剤シミュレーション．
        drugs  : createDrugDataで作成した薬剤データのリスト
        dose   : doseのリスト
        target : 標的を変える場合，リストとして標的を加える
    """
    df = pd.DataFrame()
    for index, dose in enumerate(doses):
        print("    step: {} >> ".format(index))
        drugs[0]["dose"] = dose[0]
        drugs[1]["dose"] = dose[1]
        if target:
            drugs[0]["target"] = target[0]
            drugs[1]["target"] = target[1]
        data = sim(drugs)
        data[1]["growthRate"] = [data[0] / Lambda_0]
        df = pd.concat([df, data[1]])
    df = df.reset_index(drop=True)
    drugData = pd.DataFrame([[d[0], d[1]] for d in doses], columns=["dose1", "dose2"])
    df = pd.concat([drugData, df], axis=1)
    return df

def sim_oldeval(drugName, doses, Lambda_0, target=None):
    """
        論文の評価関数を使ったシミュレーション．
        drugName : 薬剤の名前のリスト
        doses    :
        target   : 
    """
    drug1 = [createDrugData(drugName[0])]
    drug2 = [createDrugData(drugName[1])]
    drug3 = [createDrugData(drugName[0]), createDrugData(drugName[1])]
    if target:
        drug1[0]["target"] = target[0]
        drug2[0]["target"] = target[1]
        drug3[0]["target"] = target[0]
        drug3[1]["target"] = target[1]
    resultList = []
    for index, dose in enumerate(doses):
        print("    step: {} >> ".format(index))
        x = sim(drug1, [dose[0]])[0] / Lambda_0
        y = sim(drug2, [dose[1]])[0] / Lambda_0
        xy = sim(drug3, dose)[0] / Lambda_0
        if 0 in dose: Eps = 0
        else: Eps = checkEpsilon(x, y, xy)
        resultList.append([dose[0], dose[1], x, y, xy, Eps])
    data = pd.DataFrame(resultList, columns=["dose1", "dose2", "growthRate1", "growthRate2", "growthRate3", "epsilon"])
    return data

def sim_neweval(drug, slope, doseList, Lambda_0, target=None):
    if target:
        for i in range(len(target)):
            drug[i]["target"] = target[i]
    resultList = []
    for index, doses in enumerate(doseList):
        growthList = [sim(drugs, dose)[0] / Lambda_0 for dose in doses]
        resultList.append([(index + 1) * 10, slope, growthList, checkLinerType(growthList)])
    data = pd.DataFrame(resultList, columns=["MidPoint", "Slope", "resultList", "LinerType"])
    return data




if __name__ == "__main__":
    # 薬剤リスト
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]

    # IC30の計算
    makedir("IC30")
    IC30_file = "IC30/ribo8.csv"
    try:
        IC30_df = pd.read_csv(IC30_file)
        IC30 = {i: IC30_df[i][0] for i in IC30_df.columns}
    except:
        IC30 = calcIC(drugNames, {drugName: 20 for drugName in drugNames}, .3)
        IC30_df = pd.DataFrame({i: [IC30[i]] for i in drugNames})
        IC30_df.to_csv(IC30_file, index=False)
    
    Lambda_0 = sim([], [])[0]

    ## 単剤のシミュレーション
    # csvdir = "results/ribo8/csv/single"
    # makedir(csvdir)
    # result = []
    # print("start simulation >>")
    # num = int(sys.argv[-1])
    # drugName = drugNames[num]
    # print("{} >> ".format(drugName))
    # drugs = [createDrugData(drugName)]
    # doses = np.linspace(0, IC30[drugName] * 2, 101)
    # df = pd.DataFrame()
    # for index, dose in enumerate(doses):
    #     print("    step: {} >> ".format(index))
    #     drugs[0]["dose"] = dose
    #     data = sim(drugs)
    #     data[1]["growthRate"] = [data[0] / Lambda_0]
    #     df = pd.concat([df, data[1]])
    # df = df.reset_index(drop=True)
    # drugData = pd.DataFrame([[d] for d in doses], columns=["dose"])
    # df = pd.concat([drugData, df], axis=1)
    # df.to_csv("{}/{}.csv".format(csvdir, drugName), index=False)

    ## 組合せシミュレーション
    # csvdir = "results/ribo8/csv/double/normal"
    # makedir(csvdir)
    # drugNameList = itr.combinations_with_replacement(drugNames, 2)
    # num = int(sys.argv[-1])
    # print("start combination >> ")
    # for drugName in drugNameList:
    #     dirName = "{}/{}".format(csvdir, "_".join(drugName))
    #     print("{} vs {}".format(drugName[0], drugName[1]))
    #     drugs = [createDrugData(drugName[0]), createDrugData(drugName[1])]
    #     doses = divideDoses(drugName, IC30, num, 101, 101)
    #     df = sim_comb(drugs, doses, Lambda_0)
    #     df.to_csv("{}/{}.csv".format(dirName, num), index=False)

    ## 組合せシミュレーション(仮想薬剤)
    # csvdir = "results/ribo8/csv/double/virtual"
    # makedir(csvdir)
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["A", "A"], ["A", "B"], ["A", "C"]]
    # num = int(sys.argv[-1])
    # print("start combination >> ")
    # for drugName in drugNameList:
    #     print("  {} vs {} >> ".format(drugName[0], drugName[1]))
    #     doses = divideDoses(drugName, IC30, num, 101, 101)
    #     for target in targetList:
    #         dirName = "{}/{}".format(csvdir, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
    #         makedir(dirName)
    #         print("    {} vs {} >> ".format(target[0], target[1]))
    #         drugs = [createDrugData(drugName[i]) for i in range(len(drugName))]
    #         df = sim_comb(drugs, doses, Lambda_0, target)
    #         df.to_csv("{}/{}.csv".format(dirName, num), index = False)

    ## neweval simulation
    csvdir = "results/ribo8/csv/neweval/normal"
    makedir(csvdir)
    num = int(sys.argv[-1])
    drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    slopeList = [1./4., 1./2., 1., 2., 4.]
    print("start simulation >> ")
    for drugName in drugNameList:
        dirName = "{}/{}".format(csvdir, "_".join(drugName))
        makedir(dirName)
        print("  {} vs {} >>".format(drugName[0], drugName[1]))
        doses = divideDosesNeweval(drugName, IC30, num, 5, 11)
        drugs = [createDrugData(drugName[0]), createDrugData(drugName[1])]
        df = sim_neweval(drugs, slopeList[num], doses, Lambda_0)
        df.to_csv("{}/{}.csv".format(dirName, num), index=False)

    ## neweval simulation (virtual drug)
    # csvdir = "results/ribo8/csv/neweval/virtual"
    # makedir(csvdir)
    # num = int(sys.argv[-1])
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["A", "A"], ["A", "B"], ["A", "C"]]
    # slopeList = [1./4., 1./2., 1., 2., 4.]
    # print("start combination >> ")
    # for drugName in drugNameList:
    #     print("  {} vs {} >> ".format(drugName[0], drugName[1]))
    #     doses = divideDosesNeweval(drugName, IC30, num)
    #     drugs = [createDrugData(drugName[0]), createDrugData(drugName[1])]
    #     for target in targetList:
    #         print("    {} vs {} >> ".format(target[0], target[1]))
    #         dirName = "{}/{}".format(csvdir, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
    #         makedir(dirName)
    #         df = sim_neweval(drugs, slopeList[num], doses, Lambda_0, target)
    #         df.to_csv("{}/{}.csv".format(dirName, num), index=False)
