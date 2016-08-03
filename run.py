#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True

@reaction_rules
def r30_binding_reaction(a_ex, a, r30_b, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a + r30_u > r30_b | K_on * a * (r_u - r_min)
    r30_b > a + r30_u | K_off * r30_b
    a > ~a | a * Lambda # diffusion
    r30_b > ~r30_b | r30_b * Lambda # diffusion


@reaction_rules
def r50_binding_reaction(a_ex, a, r50_b, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a + r50_u > r50_b | K_on * a * (r_u - r_min)
    r50_b > a + r50_u | K_off * r50_b # dissociation
    a > ~a | a * Lambda # diffusion
    r50_b > ~r30_b | r50_b * Lambda # diffusion


@reaction_rules
def ribo_binding_reaction(a_ex, a, r_b, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a + r_u > r_b | K_on * a * (r_u - r_min)
    r_b > a + r_u | K_off * r_b
    r_b > a + r30_u + r50_u | Kd * r_b # dissociation
    a > ~a | a * Lambda # diffusion
    r_b > ~r_b | r_b * Lambda # diffusion


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
    薬剤固有の値はIC50_a, IC50, Lambda_0_a
    複数入ってきた時の対応はどうするか
    """

    Delta_r = r_max - r_min # µM
    K_off = K_on * K_D # riboと薬剤との結合
    r_u_0 = Lambda_0 / K_t + r_min
    Ka = (Lambda_0 + Kd) / (p * p * r_u_0)


    with reaction_rules():
        ### expression
        Lambda = (r_u - r_min) * K_t
        SUP = (Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))) * (1 + p) # subunit production expression

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
        r_u > r30_u + r50_u | Kd * r_u

        # diffusion
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
    # 0 : 以前のモデルとの比較
    # 1 :
    # 2

    # 保存用ディレクトリの作成
    savedir = "./images/result4"
    makedir(savedir)


    # 初期値
    r_min = 19.3 # Lambdaの計算で使用
    K_t = 6.1 * 10 ** -2 # Lambdaの計算で使用
    Kd = 6. # r30とr50の解離率
    p = .1 # r30とr50の比率

    ## drug data
    dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    Lambda_0 =  [1.35, 0.85, 0.40] # 1/h (1.35, 0.85, 0.40)
    Lambda_0_a = {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83} # 1/h
    IC50 = {"Streptmycin": [0.41, 0.28, 0.196], "Kanamycin": [0.246, 0.096, 0.065], "Tetracycline": [0.5, 0.6, 1.45], "Chloramphenicol": [2.85, 2.65, 5.7]} # µg/ml
    IC50_a = {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49} # µg/ml
    a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}
    dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}

    # 1. 以前のモデルと今回のモデルの比較
    """
    ## ribo3
    ### Kanamycinで比較
    name = "Kanamycin"
    dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}
    legend = ["r_u"]
    savename = savedir + "/ribo_ribo3.png"

    drugs = [makeDrugDatas(name, 0)]
    drug_types = ["30s", "50s", "ribo"]
    single_result = [] # 単剤用の結果の格納場所
    for num, drug_type in enumerate(drug_types):
        drugs[0]["type"] = drug_type

        # doseを振り、モデルを実行
        for count, dose in enumerate(np.linspace(0, a_ex[name], 51)):
            drugs[0]["dose"] = dose
            result, legend = run(drugs, step=100, inpData=dataset, legend=legend)
            result = (result[-1][1] - r_min) * K_t / Lambda_0[0] # 結果をgrowthに書き換え
            if num == 0:
                single_result.append([dose, result])
            else:
                single_result[count].append(result)

    ## ribo.py
    from ribo import run as run1
    dataset = {"Lambda_0": Lambda_0[0], "Lambda_0_a": Lambda_0_a[name], "IC50": IC50[name][0], "IC50_a":IC50_a[name]}
    for index, dose in enumerate(np.linspace(0, a_ex[name], 51)):
        result = run1(dose, dataset=dataset)
        single_result[index].append(result["growth"] / Lambda_0[0])

    ## make Graph
    legend = ["30s", "50s", "ribo", "previous"]
    xlabel = "Extracellular Antibiotic Concentration $a_{ex}$ ($\mu$M)"
    ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"
    title = "Single Drug Reaction"
    makeGraph(np.array(single_result), savename, legend, title, xlabel, ylabel)
    """

    # 2. nature2006のモデルとの比較
    """
    savename = savedir + "/double_bar.png"
    dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}

    ## Chloramphenicol 分子量: 323.132 g/mol, Dose: 1 µg/ml
    Chloramphenicol = 1. * 10 ** 3 / 323.132 # µM
    ## Tetracycline 分子量: 444.435 g/mol, Dose: 2 µg/ml
    Tetracycline = 2. * 10 ** 3 / 444.435 #µM
    drug_data = {"Streptmycin": {"dose": 5., "type": "30s"},
                 "Chloramphenicol": {"dose": Chloramphenicol, "type": "50s"},
                 "Tetracycline": {"dose": Tetracycline, "type": "30s"}
                }

    drug_comb = [["Chloramphenicol", "Tetracycline"],
                 ["Chloramphenicol", "Streptmycin"],
                 ["Tetracycline", "Streptmycin"]
                ]
    legend = ["r_u"]

    for index, drug in enumerate(drug_comb):
        # φ
        second_result = [1.0]
        name1 = drug[0]
        name2 = drug[1]

        # 1つめ
        drugs = [{"name": name1, "type": drug_data[name1]["type"], "dose": drug_data[name1]["dose"], "Lambda_0_a": Lambda_0_a[name1], "IC50": IC50[name1][0], "IC50_a": IC50_a[name1]}]
        result, legend = run(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
        second_result.append(result)

        # 2つめ
        drugs = [{"name": name2, "type": drug_data[name2]["type"], "dose": drug_data[name2]["dose"], "Lambda_0_a": Lambda_0_a[name2], "IC50": IC50[name2][0], "IC50_a": IC50_a[name2]}]
        result, legend = run(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
        second_result.append(result)

        # 組み合わせ
        drugs = [{"name": name1, "type": drug_data[name1]["type"], "dose": drug_data[name1]["dose"], "Lambda_0_a": Lambda_0_a[name1], "IC50": IC50[name1][0], "IC50_a": IC50_a[name1]},
                 {"name": name2, "type": drug_data[name2]["type"], "dose": drug_data[name2]["dose"], "Lambda_0_a": Lambda_0_a[name2], "IC50": IC50[name2][0], "IC50_a": IC50_a[name2]}
                ]
        result, legend = run(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
        second_result.append(result)

        # グラフ化
        x = range(4)
        plt.subplot(330 + index + 4)
        plt.bar(x, second_result, align="center")
        plt.xticks(x, ["$\phi$", name1[:3], name2[:3], "%s+\n%s" % (name1[:3], name2[:3])])
        plt.xlabel("Drug Name")
        plt.ylabel("Normalized Growth Rate $\lambda/\lambda_{0}$")
        plt.ylim(0, 1.1)
        plt.title("%s & %s" % (name1[:3], name2[:3]))
    plt.tight_layout()
    plt.savefig(savename, bbox_inches="tight", pad_inches=0.3, dpi=200)
    plt.close()
    """

    # 3. 追加データ：4つの薬剤のヒートマップ
    ##　IC50の取得
    """
    from collections import defaultdict
    dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}
    legend = ["r_u"]

    newIC50 = {} # {drugName: [dose, result], ...}
    single_result = defaultdict(list) # {drugName: [[dose, result], [dose, result] ... ]}
    for index, dName in enumerate(dNames):
        drugs = [makeDrugDatas(dName, 0)]
        # doseを振り、モデルを実行
        print dName
        for dose in np.linspace(.0, a_ex[dName], 201): # 本来はここをそれぞれのa_exにする。
            drugs[0]["dose"] = dose
            result, legend = run(drugs, step=100, inpData=dataset, legend=legend)
            result = (result[-1][1] - r_min) * K_t / Lambda_0[0] # 結果をgrowthに書き換え
            if not newIC50.get(dName):
                newIC50[dName] = [dose, result]
            else:
                if abs(0.5 - newIC50[dName][1]) > abs(0.5 - result): # IC50をnomalize
                    newIC50[dName] = [dose, result]
            single_result[dName].append([dose, result])

    ## データの整形
    ### result data
    data = []
    for index, dName in enumerate(dNames):
        for i, (x, y) in enumerate(single_result[dName]):
            if index == 0:
                data.append([x, y])
            else:
                data[i].append(x)
                data[i].append(y)

    ### newIC50
    newIC50 = {key: value[0] for key, value in newIC50.items()} # {drugName: dose, ...}
    print newIC50

    ## makeGraph
    xlabel = "Extracellular Antibiotic Concentration $a_{ex}$ ($\mu$M)"
    ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"
    title = "Single Drug Reaction"
    outputType = "multi"
    makeGraph(np.array(data), savename, dNames, title, xlabel, ylabel, outputType)
    """

    ## ヒートマップの作成
    """
    import seaborn as sns
    import pandas as pd
    import itertools as itr

    newIC50 = {'Kanamycin': 0.27250000000000002, 'Streptmycin': 0.58799999999999997, 'Chloramphenicol': 3.8000000000000003, 'Tetracycline': 0.70000000000000007}

    dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}
    legend = ["r_u"]
    savename = savedir + "/drugheatmap.png"

    # for drug in dNames:
    #     print drug
    #     drugs = [makeDrugDatas(drug, 0)]
    #     drugs[0]["dose"] = newIC50[drug]
    #     result, legend = run(drugs, inpData=dataset, legend=legend)
    #     result = (result[-1][1] - r_min) * K_t / Lambda_0[0] # 結果をgrowthに書き換え
    #     print result


    drugList = list(itr.combinations(dNames, 2)) # drugListの組み合わせを作成
    drugList = [["Kanamycin", "Streptmycin"], ["Streptmycin", "Kanamycin"]]

    # drugs = [makeDrugDatas(drugList[0][0], 0), makeDrugDatas(drugList[0][1], 0)]
    drugs = [makeDrugDatas(drugList[0][0], 0)]

    drugs[0]["dose"] = newIC50[drugList[0][0]] * 1.5
    # drugs[1]["dose"] = 0.
    result, legend = run(drugs, inpData=dataset, legend=legend)
    result = (result[-1][1] - r_min) * K_t / Lambda_0[0] # 結果をgrowthに書き換え

    print result

    # for index, drug in enumerate(drugList):
    #     drugs = [makeDrugDatas(drug[0], 0), makeDrugDatas(drug[1], 0)]
    #
    #     X = np.linspace(0, newIC50[drug[0]], 5) # 1個目の薬剤のdose
    #     Y = np.linspace(0, newIC50[drug[1]], 5) # 2個目の薬剤のdose
    #     doses = [[round(x, 4), round(y, 4)] for x in X for y in Y] # dose listの作成
    #     X = [x[0] for x in doses]
    #     Y = [y[1] for y in doses]
    #     heatmap_result = []
    #
    #     for dose in doses:
    #         print dose
    #         drugs[0]["dose"] = dose[0]
    #         drugs[1]["dose"] = dose[1]
    #         result, legend = run(drugs, inpData=dataset, legend=legend)
    #         result = (result[-1][1] - r_min) * K_t / Lambda_0[0] # 結果をgrowthに書き換え
    #         heatmap_result.append(result)
    #
    #     data = pd.DataFrame([X, Y, heatmap_result]).T
    #     data.columns = [drugs[0]["name"], drugs[1]["name"], "growth"]
    #     print data
    #     data_pivot = pd.pivot_table(data=data, values="growth", columns=drugs[1]["name"], index=drugs[0]["name"], aggfunc=np.mean)
    #     pltnum = 230 + index + 1
    #     plt.subplot(pltnum)
    #     sns.heatmap(data_pivot)
    #     plt.xlabel(drug[0])
    #     plt.ylabel(drug[1])
    # plt.savefig(savename, dpi=200)
    # plt.close()
    """

    # 4. テスト
    dName = "Kanamycin"
    dataset = {"K_D": .1, "Lambda_0": Lambda_0[0], "p": .1, "K_on": 15.}

    drugs = [makeDrugDatas(dName, 0)]

    data = []
    for dose in np.linspace(0, a_ex[dName], 201):
        drugs[0]["dose"] = dose
        drugs[0]["type"] = "ribo"
        result, legend = run(drugs, inpData=dataset, legend=["r_u"])
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0] # growthに変換
        data.append([dose, result])

    ## make Graph
    xlabel = "Extracellular Antibiotic Concentration $a_{ex}$ ($\mu$M)"
    ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"
    title = "Single Drug Reaction"
    savename = "%s/Kanamycin_ribo_Kon_20_KD_01.png" % (savedir)
    makeGraph(np.array(data), savename, title=title, xlabel=xlabel, ylabel=ylabel)
