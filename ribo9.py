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

# main
@reaction_rules
def transport(drug_ex, drug, P_in, P_out):
    # 薬剤の取り込み
    ~drug_ex == drug | (P_in * drug_ex, P_out * drug)

@reaction_rules
def dilution(mol, Lambda):
    # 希釈
    mol > ~mol | Lambda * mol

@reaction_rules
def bondingR2R(riboA, riboB, riboA_riboB, sub_k_a, sub_k_d, r30_tot, r_min):
    # リボソームの結合
    riboA + riboB == riboA_riboB | (sub_k_a, sub_k_d * (r30_tot - r_min) * (riboA_riboB / r30_tot))

@reaction_rules
def bondingD2R(drug, ribo, ribo_drug, K_on, K_off):
    # 薬剤とリボソームの結合
    drug + ribo == ribo_drug | (K_on, K_off)

def createDrugData(drugName):
    """
        薬剤のプロファイルを作成する関数
        すべて，Lambda_0 = 1.35

        name        : 薬剤の名前
        target      : 標的サブユニット(r30 or r50)
        type        : 薬剤の形式(A:遊離のみ, B:複合体のみ, C:両方(未実装))
        dose        : 投与量
     """

    datas = {
        "Lambda_0_a" : {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83},
        "IC50_a"     : {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49},
        "target"     : {"Streptmycin": "r30", "Kanamycin": "r30", "Tetracycline": "r30", "Chloramphenicol": "r50"},
        "type"       : {"Streptmycin": "A", "kanamycin": "A", "Tetracycline": "A", "Chloramphenicol": "B"}
    }

    drugData = {
        "name"       : drugName,
        "target"     : datas["target"][drugName],
        "type"       : datas["type"][drugName],
        "dose"       : .0,
        "Lambda_0_a" : datas["Lambda_0_a"][drugName],
        "IC50_a"     : datas["IC50_a"][drugName]
    }
    return drugData

def switchTarget(target):
    """
    target = r30 or r50
    """
    if target == "r30":
        val = "r50"
    else :
        val = "r30"
    return(val)

def createModel(drugs, K_D=1., K_on=3., sub_k_d=1., sub_p=1., Lambda_0=1.35, pattern=0):
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
    pattern : 競合するかしないか 0: 競合しない 1: 競合する
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
        bondingDrugList = [] # 薬剤がリボソームに結合する経路のリスト
        bondingRiboList = [] # リボソームのサブユニット同士の結合経路のリスト
        dilutionList = [] # 希釈のリスト
        r30_tot = (r30 + r30_r50) # 登場するr30
        for index, drug in enumerate(drugs):
            drug_ex = _eval("drug{}_ex".format(index))
            mol = _eval("{}{}".format(drug["type"], index))
            P_in = delta_r * drug["Lambda_0_a"] / 2.0 / drug["IC50_a"] # 薬剤の流入速度
            P_out = (drug["Lambda_0_a"] / 2) ** 2.0 / K_t / K_D # 薬剤の流出速度
            transport(drug_ex, mol, P_in, P_out) # 薬剤の流入経路の追加
            dilutionList.append(mol)
            if drug["type"] == "A":
                # 遊離サブユニットに結合 A + r > rA
                bondingDrugList.append([
                    mol,
                    _eval(drug["target"]),
                    _eval("{}A{}".format(drug["target"], index))
                ])
                dilutionList.append(_eval("{}A{}".format(drug["target"], index)))
                if drug["target"] == "r30":
                    # 標的が30sサブユニットのときだけ，作成したものをr30_totに追加
                    r30_tot += _eval("r30A{}".format(index))
                if len(drugs) > 1:
                    if drugs[(index + 1) % len(drugs)]["target"] == drug["target"] and drugs[(index + 1) % len(drugs)]["type"] == "A":
                        # 2剤がどちらも同じ標的の遊離リボソームサブユニットに結合する型だったとき
                        bondingDrugList.append([
                            mol,
                            _eval("{}A{}".format(drug["target"], (index + 1) % len(drugs))),
                            _eval("{}A0A1".format(drug["target"]))
                        ])
                        if index == 0:
                            dilutionList.append(_eval("{}A0A1".format(drug["target"])))
                            if drug["target"] == "r30":
                                # 標的が30sサブユニットのときだけ，作成したものをr30_totに追加
                                r30_tot += _eval("r30A0A1")

            elif drug["type"] == "B":
                # 複合体リボソームに結合 B + r30_r50 = r30B_r50
                if drug["target"] == "r30":
                    # 標的がr30の時
                    bondingDrugList.append([
                        mol,
                        _eval("r30_r50"),
                        _eval("r30B{}_r50".format(index))
                    ])
                    dilutionList.append(_eval("r30B{}_r50".format(index)))
                    r30_tot += _eval("r30B{}_r50".format(index))
                    if len(drugs) > 1:
                        if drugs[(index + 1) % len(drugs)]["type"] == drug["type"]:
                            if drugs[(index + 1) % len(drugs)]["target"] == drug["target"]:
                                # 標的が同じ
                                bondingDrugList.append([
                                    mol,
                                    _eval("r30B{}_r50".format((index + 1) % len(drugs))),
                                    _eval("r30B0B1_r50")
                                ])
                                if index == 0:
                                    r30_tot += _eval("r30B0B1_r50")
                                    dilutionList.append(_eval("r30B0B1_r50"))
                            else :
                                # 標的が異なる
                                bondingDrugList.append([
                                    mol,
                                    _eval("r30_r50B{}".format((index + 1) % len(drugs))),
                                    _eval("r30B{}_r50B{}".format(index, (index + 1) % len(drugs)))
                                ])
                                if index == 0:
                                    r30_tot += _eval("r30B{}_r50B{}".format(index, (index + 1) % len(drugs)))
                                    dilutionList.append(_eval("r30B{}_r50B{}".format(index, (index + 1) % len(drugs))))
                else:
                    bondingDrugList.append([
                        mol,
                        _eval("r30_r50"),
                        _eval("r30_r50B{}".format(index))
                    ])
                    dilutionList.append(_eval("r30_r50B{}".format(index)))
                    r30_tot += _eval("r30_r50B{}".format(index))
                    if len(drugs) > 1:
                        if drugs[(index + 1) % len(drugs)]["type"] == drug["type"]:
                            if drugs[(index + 1) % len(drugs)]["target"] == drug["target"]:
                                # 標的が同じ
                                bondingDrugList.append([
                                    mol,
                                    _eval("r30_r50B{}".format((index + 1) % len(drugs))),
                                    _eval("r30_r50B0B1")
                                ])
                                if index == 0:
                                    r30_tot += _eval("r30_r50B0B1")
                                    dilutionList.append(_eval("r30_r50B0B1"))
                            else :
                                # 標的が異なる
                                bondingDrugList.append([
                                    mol,
                                    _eval("r30B{}_r50".format((index + 1) % len(drugs))),
                                    _eval("r30B{}_r50B{}".format((index + 1) % len(drugs)), index)
                                ])
                                if index == 0:
                                    r30_tot += _eval("r30B{}_r50B{}".format((index + 1) % len(drugs)), index)
                                    dilutionList.append(_eval("r30B{}_r50B{}".format((index + 1) % len(drugs)), index))

        print(r30_tot)
        # r30_totを使って計算して，経路を作成
        Lambda = (r30_tot - r_min) * K_t * r30_r50 / r30_tot
        SUP = Lambda * r30_tot

        # リボソームサブユニットの生成
        ~r30 > r30 | SUP
        ~r50 > r50 | SUP

        # リボソームサブユニットの結合
        bondingRiboList.append([r30, r50, r30_r50])
        for riboReaction in bondingRiboList:
            bondingR2R(riboReaction[0], riboReaction[1], riboReaction[2], sub_k_a, sub_k_d, r30_tot, r_min)

        # 希釈
        dilutionList.extend([r30, r50, r30_r50])
        for dilutionMol in dilutionList:
            print(dilutionMol)
            dilution(dilutionMol, Lambda)

        # 薬剤のリボソームへの結合
        for drugReaction in bondingDrugList:
            bondingD2R(drugReaction[0], drugReaction[1], drugReaction[2], K_on, K_off)


    return get_model()


if __name__ == "__main__":
    drugs = [createDrugData("Streptmycin"), createDrugData("Streptmycin")]
    model = createModel(drugs=drugs)
    for i, rr in enumerate(model.reaction_rules()):
        print("{:03d}: {}".format(i, rr.as_string()))
