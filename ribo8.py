# coding: utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
import math
import itertools as itr
import seaborn as sns
import pandas as pd
util.decorator.SEAMLESS_RATELAW_SUPPORT = True

@reaction_rules
def createVariable(drug_ex, drug, P_in, P_out):
    # 薬剤の流入式を入れる関数
    ~drug_ex == drug | (P_in * drug_ex, P_out * drug)

def createDrugData(drugName):
    # 薬剤のプロファイルを作成する関数
    # すべて，Lambda_0 = 1.35
    datas = {
        "Lambda_0_a": {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83},
        "IC50_a": {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49},
        "target": {"Streptmycin": "A", "Kanamycin": "A", "Tetracycline": "B", "Chloramphenicol": "C"}
    }

    drugData = {
        "name": drugName,
        "target": datas["target"][drugName],
        "dose": .0,
        "Lambda_0_a": datas["Lambda_0_a"][drugName],
        "IC50_a": datas["IC50_a"][drugName]
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
        r30_tot = (r30 + r30_r50 + r30_r50C + r30_r50D + r30_r50CD +
            r30A + r30A_r50 + r30A_r50C + r30A_r50D + r30A_r50CD +
            r30B + r30B_r50 + r30B_r50C + r30B_r50D + r30B_r50CD +
            r30AB + r30AB_r50 + r30AB_r50C + r30AB_r50D + r30AB_r50CD)
        Lambda = (r30_tot - r_min) * K_t * r30_r50 / r30_tot
        SUP = Lambda * r30_tot
        for index, drug in enumerate(drugs):
            # 薬剤の流入の式を追加
            drug_ex = _eval("drug{}_ex".format(index))
            target = _eval(drug["target"])
            P_in = delta_r * drug["Lambda_0_a"] / 2.0 / drug["IC50_a"] # 薬剤の流入
            P_out = (drug["Lambda_0_a"] / 2) ** 2.0 / K_t / K_D # 薬剤の流出
            createVariable(drug_ex, target, P_in, P_out)

        # DrugA(30s)
        # 薬剤流入
        # ~A_ex == A | (P_in * A_ex, P_out)
        # 薬剤希釈
        A > ~A | Lambda * A
        # 薬剤結合
        ## 標的サブユニットへの結合
        A + r30 == r30A | (K_on, K_off)
        A + r30B == r30AB | (K_on, K_off)
        ## 機能性リボソームへの結合
        A + r30_r50 == r30A_r50 | (K_on, K_off)
        A + r30_r50C == r30A_r50C | (K_on, K_off)
        A + r30_r50D == r30A_r50D | (K_on, K_off)
        A + r30_r50CD == r30A_r50CD | (K_on, K_off)
        A + r30B_r50 == r30AB_r50 | (K_on, K_off)
        A + r30B_r50C == r30AB_r50C | (K_on, K_off)
        A + r30B_r50D == r30AB_r50D | (K_on, K_off)
        A + r30B_r50CD == r30AB_r50CD | (K_on, K_off)

        # DrugB(30s)
        # 薬剤流入
        # ~B_ex == B | (P_in * B_ex, P_out)
        # 薬剤希釈
        B > ~B | Lambda * B
        # 薬剤結合
        ## 標的サブユニットへの結合
        B + r30 == r30B | (K_on, K_off)
        B + r30A == r30AB | (K_on, K_off)
        ## 機能性リボソームへの結合
        B + r30_r50 == r30B_r50 | (K_on, K_off)
        B + r30_r50C == r30B_r50C | (K_on, K_off)
        B + r30_r50D == r30B_r50D | (K_on, K_off)
        B + r30_r50CD == r30B_r50CD | (K_on, K_off)
        B + r30A_r50 == r30AB_r50 | (K_on, K_off)
        B + r30A_r50C == r30AB_r50C | (K_on, K_off)
        B + r30A_r50D == r30AB_r50D | (K_on, K_off)
        B + r30A_r50CD == r30AB_r50CD | (K_on, K_off)

        # DrugC(50s)
        # 薬剤流入
        # ~C_ex == C | (P_in * C_ex, P_out)
        # 薬剤希釈
        C > ~C | Lambda * C
        # 薬剤結合
        ## 標的サブユニットへの結合
        C + r50 == r50C | (K_on, K_off)
        C + r50D == r50CD | (K_on, K_off)
        ## 機能性リボソームへの結合
        C + r30_r50 == r30_r50C | (K_on, K_off)
        C + r30_r50D == r30_r50CD | (K_on, K_off)
        C + r30A_r50 == r30A_r50C | (K_on, K_off)
        C + r30A_r50D == r30A_r50CD | (K_on, K_off)
        C + r30B_r50 == r30B_r50C | (K_on, K_off)
        C + r30B_r50D == r30B_r50CD | (K_on, K_off)
        C + r30AB_r50 == r30AB_r50C | (K_on, K_off)
        C + r30AB_r50D == r30AB_r50CD | (K_on, K_off)

        # DrugD(50s)
        # 薬剤流入
        # ~D_ex == D | (P_in * D_ex, P_out)
        # 薬剤希釈
        D > ~D | Lambda * D
        # 薬剤結合
        ## 標的サブユニットへの結合
        D + r50 == r50D | (K_on, K_off)
        D + r50C == r50CD | (K_on, K_off)
        ## 機能性リボソームへの結合
        D + r30_r50 == r30_r50D | (K_on, K_off)
        D + r30_r50C == r30_r50CD | (K_on, K_off)
        D + r30A_r50 == r30A_r50D | (K_on, K_off)
        D + r30A_r50C == r30A_r50CD | (K_on, K_off)
        D + r30B_r50 == r30B_r50D | (K_on, K_off)
        D + r30B_r50C == r30B_r50CD | (K_on, K_off)
        D + r30AB_r50 == r30AB_r50D | (K_on, K_off)
        D + r30AB_r50C == r30AB_r50CD | (K_on, K_off)

        # リボソームの結合パターン
        r30 + r50 == r30_r50 | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30_r50 / r30_tot) # もともとはこれだけ
        r30 + r50C == r30_r50C | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30_r50C / r30_tot)
        r30 + r50D == r30_r50D | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30_r50D / r30_tot)
        r30 + r50CD == r30_r50CD | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30_r50CD / r30_tot)
        r30A + r50 == r30A_r50 | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30A_r50 / r30_tot)
        r30A + r50C == r30A_r50C | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30A_r50C / r30_tot)
        r30A + r50D == r30A_r50D | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30A_r50D / r30_tot)
        r30A + r50CD == r30A_r50CD | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30A_r50CD / r30_tot)
        r30B + r50 == r30B_r50 | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30B_r50 / r30_tot)
        r30B + r50C == r30B_r50C | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30B_r50C / r30_tot)
        r30B + r50D == r30B_r50D | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30B_r50D / r30_tot)
        r30B + r50CD == r30B_r50CD | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30B_r50CD / r30_tot)
        r30AB + r50 == r30AB_r50 | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30AB_r50 / r30_tot)
        r30AB + r50C == r30AB_r50C | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30AB_r50C / r30_tot)
        r30AB + r50D == r30AB_r50D | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30AB_r50D / r30_tot)
        r30AB + r50CD == r30AB_r50CD | (sub_k_a, sub_k_d * (r30_tot - r_min) * r30AB_r50CD / r30_tot)

        # リボソームサブユニットの希釈
        r30 > ~r30 | Lambda * r30
        r50 > ~r50 | Lambda * r50
        r30A > ~r30A | Lambda * r30A
        r30B > ~r30B | Lambda * r30B
        r30AB > ~r30AB | Lambda * r30AB
        r50C > ~r50C | Lambda * r50C
        r50D > ~r50D | Lambda * r50D
        r50CD > ~r50CD | Lambda * r50CD
        r30_r50 > ~r30_r50 | Lambda * r30_r50
        r30_r50C > ~r30_r50C | Lambda * r30_r50C
        r30_r50D > ~r30_r50D | Lambda * r30_r50D
        r30_r50CD > ~r30_r50CD | Lambda * r30_r50CD
        r30A_r50 > ~r30A_r50 | Lambda * r30A_r50
        r30A_r50C > ~r30A_r50C | Lambda * r30A_r50C
        r30A_r50D > ~r30A_r50D | Lambda * r30A_r50D
        r30A_r50CD > ~r30A_r50CD | Lambda * r30A_r50CD
        r30B_r50 > ~r30B_r50 | Lambda * r30B_r50
        r30B_r50C > ~r30B_r50C | Lambda * r30B_r50C
        r30B_r50D > ~r30B_r50D | Lambda * r30B_r50D
        r30B_r50CD > ~r30B_r50CD | Lambda * r30B_r50CD
        r30AB_r50 > ~r30AB_r50 | Lambda * r30AB_r50
        r30AB_r50C > ~r30AB_r50C | Lambda * r30AB_r50C
        r30AB_r50D > ~r30AB_r50D | Lambda * r30AB_r50D
        r30AB_r50CD > ~r30AB_r50CD | Lambda * r30AB_r50CD

        # リボソームサブユニットの生成
        ~r30 > r30 | SUP
        ~r50 > r50 | SUP

    return get_model()

def run(drugs=[], step=50., inpData={}, y0={"r30": 30., "r50": 30., "r30_r50": 30.}, sp_list=["r30_r50"]):
    # シミュレーションを実行する関数
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
    # 返すデータはフルデータをリターンする(species_listは入力しなかった場合はどうやって出す？？)
    return data

def calcGrowthRate(r30_r50, r30_tot, Lambda_0=1.35):
    r_min = 19.3
    K_t = 6.1 * 10 ** -2
    Lambda = (r30_tot - r_min) * K_t * r30_r50 / r30_tot
    growth = Lambda / Lambda_0
    return result

def sim(drugs, step=0.5, inpData={}, y0={"r30": 30., "r50": 30., "r30_r50": 30.}, Lambda_0=1.35):
    # rtot を計算するために sp_listを作成
    sp_list = ["r30", "r30_r50", "r30_r50C", "r30_r50D", "r30_r50CD",
        "r30A", "r30A_r50", "r30A_r50C", "r30A_r50D", "r30A_r50CD",
        "r30B", "r30B_r50", "r30B_r50C", "r30B_r50D", "r30B_r50CD",
        "r30AB", "r30AB_r50", "r30AB_r50C", "r30AB_r50D", "r30AB_r50CD"]
    result = run(drugs, step, inpData, y0, sp_list=sp_list)[-1][1:] # tを除去したリストを作成
    r30_r50 = result[1]
    r30_tot = sum(result)
    growth = calcGrowthRate(r30_r50, r30_tot)
    df = pd.DataFrame([result], columns=sp_list)
    return (growth, df)

if __name__ == "__main__":
    # 薬剤リスト
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    # species list(rtotを計算するために必要)

    drugs = [createDrugData("Streptmycin")]
    sp_list = ["r30", "r30_r50", "r30_r50C", "r30_r50D", "r30_r50CD",
        "r30A", "r30A_r50", "r30A_r50C", "r30A_r50D", "r30A_r50CD",
        "r30B", "r30B_r50", "r30B_r50C", "r30B_r50D", "r30B_r50CD",
        "r30AB", "r30AB_r50", "r30AB_r50C", "r30AB_r50D", "r30AB_r50CD"]
    doses = np.linspace(0, 5, 6)
    result = []
    drugs[0]["dose"] = 2.
    df = pd.DataFrame()
    for dose in doses:
        drugs[0]["dose"] = dose
        data = sim(drugs)
        result.append(data[0])
        df = pd.concat([df, data[1]])
    df = df.reset_index(drop=True)
    drugData = pd.DataFrame([[d] for d in doses], columns=["Streptmycin"])
    df = pd.concat([drugData, df], axis=1)
    df.to_csv("results/ribo8/test/test.csv", index=False)
    print(df)
    # print(result)
    # plt.plot(doses, result)
    # plt.show()
