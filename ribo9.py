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
    ~drug_ex == drug | (P_in * drug_ex, P_out, drug)

@reaction_rules
def dilution(mol):
    # 希釈
    mol > ~mol | Lambda * mol

@reaction_rules
def combination_ribosome(r30, r50, r30_r50):
    # リボソームの結合
    r30 + r50 == r30_r50 | (sub_k_a, sub_k_d * (r30_tot - r_min) * (r30_r50 / r30_tot))

@reaction_rules
def combination_drug(drug, ribo, ribo_drug):
    # 薬剤とリボソームの結合
    drug + ribo == ribo_drug | (K_on, K_off)

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
        r30_tot = (r30 + r30_r50)
        for index, drug in enumerate(drugs):
            drug_ex = _eval("drug{}_ex".format(index))
            mol = _eval("{}{}".format(drug["target"], index))
            P_in = delta_r * drug["Lambda_0_a"] / 2.0 / drug["IC50_a"] # 薬剤の流入
            P_out = (drug["Lambda_0_a"] / 2) ** 2.0 / K_t / K_D # 薬剤の流出
            transport(drug_ex, mol, P_in, P_out)
            dilution(mol)
            if drug["target"] == "A":
                # 遊離サブユニットに結合 A + r30 > r30A
                combination_drug(mol, _eval("r30"), _eval("r30A{}".format(index)))
                if drugs[index - 1]["target"] == "A" and pattern == 0:
                    combination_drug(mol, _eval("r30A{}".format(abs(index - 1)), _eval("r30A0A1"))
            if drug["target"] == "B":
                # 複合体リボソームに結合
                combination_drug(mol, _eval("r30_r50"), _eval("r30B{}_r50".format(index)))
                if drugs[index - 1]["target"] == "B" and pattern == 0:
                    combination_drug(mol, _eval("r30B{}_r50".format(abs(index - 1)), _eval("r30B1B2_r50"))
                elif drugs[index - 1]["target"] == "D":
                    combination_drug(mol, _eval("r30_r50D{}".format(len(drugs) - 1)), _eval("r30B{}_r50D{}".format(index, len(drugs) - index)))
