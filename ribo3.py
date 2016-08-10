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


def createModel(drugs=[], r_max=65.8, r_min=19.3, K_D=1., K_t=6.1*10**-2, K_on=3.0, Lambda_0=1.35, Kd=1., p=1.):
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
    Ka = (Kd / K_t + r_u_0) * Lambda_0 / ((p * r_u_0) ** 2) #

    with reaction_rules():
        ### expression
        Lambda = (r_u - r_min) * K_t
        SUP = (Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))) * (1 + p) # subunit product expression

        ### reaction
        ## drug
        if len(drugs) > 0:
            if drugs[0]["type"] == "30s":
                # print "drug1 targets 30s ribosomal subunit >>"
                r30_binding_reaction(a1_ex, a1, r30_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
            elif drugs[0]["type"] == "50s":
                # print "drug1 targets 50s ribosomal subunit >>"
                r50_binding_reaction(a1_ex, a1, r50_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
            elif drugs[0]["type"] == "ribo":
                # print "drug1 targets ribosome >>"
                ribo_binding_reaction(a1_ex, a1, r_1_b, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)

        if len(drugs) > 1:
            if drugs[1]["type"] == "30s":
                # print "drug2 targets 30s ribosomal subunit >>"
                r30_binding_reaction(a2_ex, a2, r30_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on, K_off, Lambda)
            elif drugs[1]["type"] == "50s":
                # print "drug2 targets 50s ribosomal subunit >>"
                r50_binding_reaction(a2_ex, a2, r50_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on, K_off, Lambda)
            elif drugs[1]["type"] == "ribo":
                # print "drug2 targets ribosome >>"
                ribo_binding_reaction(a2_ex, a2, r_2_b, drugs[1]["P_in"], drugs[1]["P_out"], K_on, K_off, Lambda)

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


def makeDrugDatas(drugName, medium):
    """
    薬剤データを作成して返す関数
    drugName : 投与する薬剤の名前
    medium : 培地の番号(0: ,
                      1: ,
                      2: ,
                     )
    """
    dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
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
    inpData:
    y0:
    """

    dataset = {"r_max": 65.8, "r_min": 19.3, "K_D": 1.0, "K_t": 6.1*10**-2, "K_on": 3.0, "Lambda_0": 1.35, "Kd": 1., "p": 1.} # 基本のデータセット
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

    # reactionの確認
    # for i, rr in enumerate(model.reaction_rules()):
    #     print("{:d}: {:s}".format(i, rr.as_string()))

    # legend
    if not legend:
        legend = y0.keys()

    runsim = run_simulation(step, solver="ode", y0=y0,
                            return_type="observer", model=model,
                            species_list=legend)
    data = runsim.data()
    return data, legend


def makeGraph(data, savename, legend=[], title="", xlabel="", ylabel="", type="single"):
    """
    グラフ作成用モジュール
    """
    import sys
    if type == "single":
        for i in range(len(data[0]) - 1):
            if len(legend) > i:
                plt.plot(data.T[0], data.T[i+1], label=legend[i])
            else:
                plt.plot(data.T[0], data.T[i+1])
    elif type == "multi":
        if len(data[0]) % 2. > 0:
            sys.stderr.write("Error: Input data length is not Even.")

        else:
            for i in range(len(data[0]) / 2):
                if len(legend) > i:
                    plt.plot(data.T[i * 2], data.T[i * 2 + 1], label=legend[i])
                else:
                    plt.plot(data.T[i * 2], data.T[i * 2 + 1])

    if title: # titleがあったらつける
        plt.title(title)

    if xlabel: # xlabelがあったらつける
        plt.xlabel(xlabel)

    if ylabel: # ylabelがあったらつける
        plt.ylabel(ylabel)

    if legend:
        plt.legend(loc="upper right")
    # plt.show()
    plt.savefig(savename, dpi=200)
    plt.close()


def makedir(dirname):
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    del(os)

if __name__ == "__main__":
    # 0 : 単剤

    # 保存用ディレクトリの作成
    savedir = "./images/result1"
    makedir(savedir)

    # 初期値
    r_min = 19.3 # Lambdaの計算で使用
    K_t = 6.1 * 10 ** -2 # Lambdaの計算で使用
    # (result[-1][1] - r_min) * K_t / Lambda_0[0] # growth_rate / Lambda_0

    ## drug data
    dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    Lambda_0 =  [1.35, 0.85, 0.40] # 1/h (1.35, 0.85, 0.40)
    Lambda_0_a = {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83} # 1/h
    IC50 = {"Streptmycin": [0.41, 0.28, 0.196], "Kanamycin": [0.246, 0.096, 0.065], "Tetracycline": [0.5, 0.6, 1.45], "Chloramphenicol": [2.85, 2.65, 5.7]} # µg/ml
    IC50_a = {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49} # µg/ml
    a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}
    dataset = {"Lambda_0": Lambda_0[0]}

    # 0 : 単剤の出力
    plt.figure(figsize=(10, 8))
    for index, dName in enumerate(dNames):
        drugs = [makeDrugDatas(dName, 0)] # 薬剤データの作成
        data = []
        for dose in np.linspace(0, a_ex[dName] * 3, 201):
            drugs[0]["dose"] = dose
            result, legend = run(drugs, step=100, legend=["r_u"])
            result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
            data.append([dose, result])

        title = "%s Reaction" % (dName)
        xlabel = "Extracellular Antibiotic Concentration $a_{ex}$ ($\mu$M)"
        ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"

        plt.subplot(220 + index + 1)
        data = np.array(data)
        plt.plot(data.T[0], data.T[1])
        plt.title(title)
        plt.axis(fontsize=8)
        plt.xlim(0, a_ex[dName] * 3)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    savename = "%s/single.png" % (savedir)
    plt.tight_layout()
    plt.savefig(savename, dpi=200)
    plt.close()

    # 1 : 2分法でIC50を計算し、ダブルノックダウンを出力
    ## 2分法でIC50を計算
    """
    newIC50 = {}
    for dName in dNames:
        drugs = [makeDrugDatas(dName, 0)] # 薬剤データの作成
        dose_max = a_ex[dName] * 3.
        dose_min = 0.
        result = 1.
        dose = 0
        while abs(0.5 - abs(result)) > 0.01:
            dose = (dose_max + dose_min) / 2. # maxとminの中間を取る
            drugs[0]["dose"] = dose
            result, legend = run(drugs, step=100, legend=["r_u"])
            result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
            if result < 0.5:
                dose_max = dose # 0.5より小さい時は
            else:
                dose_min = dose
        newIC50[dName] = dose

    newIC50 = {'Kanamycin': 0.67584228515625, 'Streptmycin': 1.458984375, 'Chloramphenicol': 10.3125, 'Tetracycline': 2.0625}

    ## ダブルノックダウン
    import itertools as itr
    drug_comb = list(itr.combinations(dNames, 2))
    data = []
    legendList = []
    for index, dList in enumerate(drug_comb):
        drugs = [makeDrugDatas(dList[0], 0), makeDrugDatas(dList[1], 0)]
        legendList.append("%s & %s" % (dList[0][:3], dList[1][:3]))
        for num, i in enumerate(np.linspace(0, 100, 101)):
            drugs[0]["dose"] = newIC50[dList[0]] * i / 100
            drugs[1]["dose"] = newIC50[dList[1]] * i / 100
            result, legend = run(drugs, step=100, legend=["r_u"])
            result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
            if index == 0:
                data.append([i, result])
            else:
                data[num].append(result)
    savename = "%s/double.png" % (savedir)
    title = "double Reaction"
    xlabel = "Ratio of Dose in Comparison with Normalized IC50 (%)"
    ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"
    makeGraph(np.array(data), savename, legendList, title, xlabel, ylabel)
    """

    # 2 : ヒートマップ作成
    """
    import seaborn as sns
    import pandas as pd
    import itertools as itr

    newIC50 = {'Kanamycin': 0.67584228515625, 'Streptmycin': 1.458984375, 'Chloramphenicol': 10.3125, 'Tetracycline': 2.0625}

    plt.figure(figsize=(10, 6))
    drug_comb = list(itr.combinations(dNames, 2)) # 薬剤の組み合わせ
    for index, dList in enumerate(drug_comb):
        drugs = [makeDrugDatas(dList[0], 0), makeDrugDatas(dList[1], 0)]
        X = np.linspace(0, newIC50[dList[0]], 5)
        Y = np.linspace(0, newIC50[dList[1]], 5)
        doses = list(itr.product(X, Y))
        for i, doseList in enumerate(doses):
            drugs[0]["dose"] = doseList[0]
            drugs[1]["dose"] = doseList[1]
            result, legend = run(drugs, step=100, legend=["r_u"])
            result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
            if i == 0:
                data = pd.DataFrame([[round(doseList[0], 3), round(doseList[1], 3), result]], columns=[dList[0], dList[1], "growth_rate"])
            else:
                data = data.append(pd.DataFrame([[doseList[0], doseList[1], result]], columns=[dList[0], dList[1], "growth_rate"]))
        heatmap = pd.pivot_table(data=data, values="growth_rate", index=dList[0], columns=dList[1], aggfunc=np.mean) # x軸を0, y軸を1番目の薬剤にしてグラフデータ化
        plt.subplot(230 + index + 1) # 1つの画像データに集約
        sns.heatmap(heatmap)
        plt.axis(fontsize=8)
        plt.xlabel(dList[0])
        plt.ylabel(dList[1])

    savename = "%s/heatmap.png" % (savedir)
    plt.fontsize()
    plt.tight_layout()
    plt.savefig(savename, dpi=200)
    plt.close()
    """

    # 4 : bar plot
