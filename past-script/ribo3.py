#coding:utf-8
# test
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

def makeCmap(color_range={"red": 30, "pink": 5, "white": 10, "light_green": 5, "green": 13, "blue": 17},
            color_num = {"red": "#ff0000", "pink": "#ffc0cb", "white": "#ffffff", "light_green": "#90ee90", "green": "#008000", "blue": "#0000ff"},
            color_list = ["red", "pink", "white", "light_green", "green", "blue"]): # eをかませるためのカラーマップを作る関数

    import matplotlib

    color_result = []
    for color in color_list:
        for i in range(color_range[color]):
            color_result.append(color_num[color])
    cmap = matplotlib.colors.ListedColormap(color_result)
    del(matplotlib)
    return cmap

def epsilon(x, y, val):
    """
    Nature genetics 2006のepsilonを計算する関数
    e = (Wxy - WxWy)/|min(Wx, Wy) - WxWy|
    """
    result = (val - x * y) / abs(min(x, y) - x * y)
    return result


if __name__ == "__main__":
    # 保存用ディレクトリの作成
    savedir = "./images/result2"
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

    # 0 : 単剤出力の比較
    """
    ## ribo.py
    import ribo

    ### IC50を取得
    oldIC50 = {}
    for dName in dNames:
        dataset = {"Lambda_0": Lambda_0[0],
                   "Lambda_0_a": Lambda_0_a[dName],
                   "IC50": IC50[dName][0],
                   "IC50_a": IC50_a[dName]}
        dose_max = a_ex[dName] * 3.
        dose_min = 0.
        result = 1.
        dose = 0
        while abs(0.5 - abs(result)) > 0.01:
            dose = (dose_max + dose_min) / 2. # maxとminの中間を取る
            result = ribo.run(dose, dataset=dataset, step=100)
            result = result["growth"] / Lambda_0[0]
            if result < 0.5:
                dose_max = dose # 0.5より小さい時は
            else:
                dose_min = dose
        oldIC50[dName] = dose

    ### 出力結果をIC50で正規化してデータ化
    data_ribo = {}
    for dName in dNames:
        result_list = []
        dataset = {"Lambda_0": Lambda_0[0],
                   "Lambda_0_a": Lambda_0_a[dName],
                   "IC50": IC50[dName][0],
                   "IC50_a": IC50_a[dName]}
        for dose in np.linspace(0, a_ex[dName], 201):
            result = ribo.run(dose, dataset=dataset, step=100)
            result = result["growth"] / Lambda_0[0]
            result_list.append([dose / oldIC50[dName], result])
        data_ribo[dName] = np.array(result_list)

    ## ribo3.py
    ### IC50を取得
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

    ### 出力結果をIC50で正規化してデータ化
    data_ribo3 = {} # 正規化したデータを薬剤名をkeyにして保存
    for index, dName in enumerate(dNames):
        result_list = []
        drugs = [makeDrugDatas(dName, 0)] # 薬剤データの作成
        for dose in np.linspace(0, a_ex[dName] * 3.5, 201):
            drugs[0]["dose"] = dose
            result, legend = run(drugs, step=100, legend=["r_u"])
            result = (result[-1][1] - r_min) * K_t / Lambda_0[0] # growth_rate
            result_list.append([dose / newIC50[dName], result]) # IC50で正規化
        data_ribo3[dName] = np.array(result_list)

    ## riboとribo3の結果をまとめてグラフ化
    plt.figure(figsize=(10, 8))
    for index, dName in enumerate(dNames):
        plt.subplot(220 + index + 1)
        plt.plot(data_ribo[dName].T[0], data_ribo[dName].T[1], label="Previous model")
        plt.plot(data_ribo3[dName].T[0], data_ribo3[dName].T[1], label="Latest model")
        plt.title("%s" % (dName))
        plt.xlabel("$a_{ex}$ conc. Normalised IC50")
        plt.ylabel("Normalized Growth Rate $\lambda/\lambda_{0}$")
        if data_ribo[dName].T[0][-1] < data_ribo3[dName].T[0][-1]:
            plt.xlim(0, data_ribo3[dName].T[0][-1])
        else:
            plt.xlim(0, data_ribo[dName].T[0][-1])
        plt.legend(loc="upper right")
    savename = "%s/single.png" % (savedir)
    plt.tight_layout()
    plt.savefig(savename, dpi=200)
    plt.close()
    """

    # 1 : ヒートマップの作成
    ## 2分法でIC50を計算
    newIC50 = {}
    for dName in dNames:
        drugs = [makeDrugDatas(dName, 0)] # 薬剤データの作成
        dose_max = a_ex[dName] * 3.
        dose_min = 0.
        result = 1.
        dose = 0
        target = 0.5
        while abs(target - abs(result)) > 0.01:
            dose = (dose_max + dose_min) / 2. # maxとminの中間を取る
            drugs[0]["dose"] = dose
            result, legend = run(drugs, step=100, legend=["r_u"])
            result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
            if result < target:
                dose_max = dose # targetより小さい時は
            else:
                dose_min = dose
        newIC50[dName] = dose

    ## ヒートマップ作成
    import seaborn as sns
    import pandas as pd
    import itertools as itr

    cmap = makeCmap() # カラーマップを設定

    # newIC50 = {'Kanamycin': 0.67584228515625, 'Streptmycin': 1.458984375, 'Chloramphenicol': 10.3125, 'Tetracycline': 2.0625}

    plt.figure(figsize=(15, 9))
    drug_comb = list(itr.combinations(dNames, 2)) # 薬剤の組み合わせ
    for index, dList in enumerate(drug_comb):
        drugs = [makeDrugDatas(dList[0], 0), makeDrugDatas(dList[1], 0)]
        X = np.linspace(0, newIC50[dList[0]]*2, 10) # IC50にするときは2倍
        Y = np.linspace(0, newIC50[dList[1]]*2, 10) # IC50にするときは2倍
        doses = list(itr.product(X, Y))
        for i, doseList in enumerate(doses):
            if (doseList[0] == X[0] or doseList[0] == X[-1]) or (doseList[1] == Y[0] or doseList[1] == Y[-1]):
                result = 0.
            else:
                drugs[0]["dose"] = doseList[0]
                drugs[1]["dose"] = doseList[1]
                result, legend = run([drugs[0]], step=100, legend=["r_u"])
                x = (result[-1][1] - r_min) * K_t / Lambda_0[0]
                result, legend = run([drugs[1]], step=100, legend=["r_u"])
                y = (result[-1][1] - r_min) * K_t / Lambda_0[0]
                result, legend = run(drugs, step=100, legend=["r_u"])
                result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
                result = epsilon(x, y, result)
            if i == 0:
                data = pd.DataFrame([[round(doseList[0], 1), round(doseList[1], 1), result]], columns=[dList[0], dList[1], "growth_rate"])
            else:
                data = data.append(pd.DataFrame([[round(doseList[0], 1), round(doseList[1], 1), result]], columns=[dList[0], dList[1], "growth_rate"]))
        heatmap = pd.pivot_table(data=data, values="growth_rate", index=dList[0], columns=dList[1], aggfunc=np.mean) # x軸を0, y軸を1番目の薬剤にしてグラフデータ化
        plt.subplot(230 + index + 1) # 1つの画像データに集約
        sns.heatmap(heatmap, vmin=-2, vmax=2, cmap=cmap, annot=True)
        # sns.heatmap(heatmap)
        plt.axis(fontsize=3)
        plt.ylabel(dList[0])
        plt.xlabel(dList[1])

    savename = "%s/epsilon_IC50.png" % (savedir)
    plt.tight_layout()
    plt.savefig(savename, dpi=200)
    plt.close()

    # 3 : bar plot
    """
    drug_comb = [["Chloramphenicol", "Tetracycline"],
                 ["Chloramphenicol", "Streptmycin"],
                 ["Tetracycline", "Streptmycin"]]
    target_growth = [{"Chloramphenicol": 0.39, "Tetracycline": 0.70},
                     {"Chloramphenicol": 0.43, "Streptmycin": 0.1}, # Streptmycin: 0.03
                     {"Tetracycline": 0.73, "Streptmycin": 0.70}
                    ]
    for index, drugList in enumerate(drug_comb):
        doses = {}
        ## 二分法をもちいて、該当するgrowth_rateを出すdoseを算出
        for dName in drugList:
            drugs = [makeDrugDatas(dName, 0)] # 薬剤データの作成
            dose_max = a_ex[dName] * 3.
            dose_min = 0.
            result = 1.
            dose = 0
            count = 0
            while abs(target_growth[index][dName] - abs(result)) > 0.001 or count < 100:
                dose = (dose_max + dose_min) / 2. # maxとminの中間を取る
                drugs[0]["dose"] = dose
                result, legend = run(drugs, step=100, legend=["r_u"])
                result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
                if result < target_growth[index][dName]:
                    dose_max = dose # 標的のgrowthより小さい時は
                else:
                    dose_min = dose
                count += 1
            doses[dName] = dose

        print doses

        ##グラフ化用データ作成
        ### φ
        data = [1.0]

        ### 1つめ
        drugs = [makeDrugDatas(drugList[0], 0)]
        drugs[0][dose] = doses[drugList[0]]
        result, legend = run(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
        data.append(result)

        ### 2つめ
        drugs = [makeDrugDatas(drugList[1], 0)]
        drugs[0][dose] = doses[drugList[1]]
        result, legend = run(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
        data.append(result)

        ### 組み合わせ
        drugs = [makeDrugDatas(drugList[0], 0),
                 makeDrugDatas(drugList[1], 0)]
        drugs[0][dose] = doses[drugList[0]]
        drugs[1][dose] = doses[drugList[1]]
        result, legend = run(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
        data.append(result)

        ## データ作成
        x = range(4)
        plt.subplot(130 + index + 1)
        plt.bar(x, data, align="center")
        plt.xticks(x, ["$\phi$", drugList[0][:3], drugList[1][:3], "%s+\n%s" % (drugList[0][:3], drugList[1][:3])])
        plt.xlabel("Drug Name")
        plt.ylabel("Normalized Growth Rate $\lambda/\lambda_{0}$")
        plt.ylim(0, 1.1)
        plt.title("%s & %s" % (drugList[0][:3], drugList[1][:3]))
    plt.tight_layout()
    savename = savedir + "/bar.png"
    plt.savefig(savename, bbox_inches="tight", pad_inches=0.3, dpi=200)
    plt.close()
    """
