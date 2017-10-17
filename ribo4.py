#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True
import math
import itertools as itr
import seaborn as sns
import pandas as pd

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


def run(drugs=[], step=50., medium=0, legend=[], inpData={}, y0={"r30_u": 30., "r50_u": 30., "r_u": 30., "r_b": .0}):
    """
    ribosomeモデル実行関数

    drugs: [{"name":, "type":, "dose":, }]
    step:
    legend:
    inpData:
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

### modules ###
def sim(drugs, dose):
    # runなどのシミュレーションをして，結果を返す関数
    for i in range(len(drugs)):
        drugs[i]["dose"] = dose[i]
    result, legend = run(drugs, legend=["r_u"])
    return calcGrowthRate(result[-1][1])

def makedir(dirname):
    """ディレクトリ作成関数"""
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    del(os)

def divideFigure(drugNames):
    """
        drugNamesに応じて，figureの分割数を計算するモジュール．
    """
    count = 0
    for i in drugNames: count += 1
    x = int(math.ceil(math.sqrt(count)))
    y = int(math.ceil(count / x))
    return (x, y)

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

def generate_cmap(colors):
    """自分で定義したカラーマップを返す(線形補完)"""
    from matplotlib.colors import LinearSegmentedColormap
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


def checkEpsilon(dName, dose): # nature genesis 2006's evaluation
    x = sim([makeDrugDatas(dName[0])], [dose[0]])
    y = sim([makeDrugDatas(dName[1])], [dose[1]])
    xy = sim([makeDrugDatas(dName[0]), makeDrugDatas(dName[1])], dose)
    return (xy - x * y) / abs(min(x, y) - x * y)

def calcGrowthRate(a, r_min=19.3, K_t=6.1*10**-2, Lambda_0=1.35):
    """
    growth_rateを計算する関数
    """
    result = (a - r_min) * K_t / Lambda_0
    return result


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
            result = calcGrowthRate(result[-1][1])
            if result < target:
                dose_max = dose
            else:
                dose_min = dose
        calc_result[dName] = dose
    return calc_result

def createHeatmap(data, drugNames, cbar=False, cmap=False, axizFontSize=16, labelSize=16):
    """
        function of create Heatmap.
        data: pandas data.
        cmap: color map of heatmap.
    """
    if not cmap: cmap = sns.diverging_palette(220, 10, as_cmap=True) # coler map
    heatmapData = pd.pivot_table(data=data, values="growth", index="a1", columns="a2") # heatmap data
    ax = sns.heatmap(heatmapData, cbar=cbar, cmap=cmap) # create Heatmap
    ax.invert_yaxis()

    a1DoseList = list(set(data["a1"].tolist()))[::2] # y
    a2DoseList = list(set(data["a2"].tolist()))[::2] # x

    ax.set_xticks(list(np.linspace(0.5, 10.5, len(a1DoseList))))
    ax.set_xticklabels(list(map(str, a2DoseList)))
    ax.set_yticks(list(np.linspace(0.5, 10.5, len(a1DoseList))))
    ax.set_yticklabels(list(map(str, a1DoseList)))

    ax.set_ylabel(drugNames[0], fontsize=axizFontSize) # create ylabel
    ax.set_xlabel(drugNames[1], fontsize=axizFontSize) # create xlabel

    ## virtual drug
    # if drugNames[0] == "Streptmycin": ax.set_ylabel("Pattern A", fontsize=axizFontSize) # create ylabel
    # else : ax.set_ylabel("Pattern B", fontsize=axizFontSize) # create ylabe
    # if drugNames    [1] == "Streptmycin": ax.set_xlabel("Pattern A", fontsize=axizFontSize) # create xlabel
    # else: ax.set_xlabel("Pattern B", fontsize=axizFontSize) # create xlabel

    ax.tick_params(labelsize=labelSize)



def oldeval(dNameList, IC30, subplot, slopeList=[1./4, 1./2, 1., 2., 4.], titleFontSize=16, axizFontSize=16, labelSize=16, csvdir="."):
    # 論文の評価関数を使用してヒートマップを作成する関数
    figsize = (subplot[1] * 10, subplot[0] * 10)
    if figsize: fig = plt.figure(figsize=figsize)
    else: fig = plt.figure()
    for index, name in enumerate(dNameList):
        midPointList = [x for x in itr.zip_longest(np.linspace(0, IC30[name[0]], 11), np.linspace(0, IC30[name[1]], 11))][1:] # 0, 0を除いた10点．中点なので，IC30でOK．
        doseList = [[[midPoint[0] * (1 + slope), midPoint[1] * (1 + (1 / slope))] for midPoint in midPointList] for slope in slopeList] # doseの組合せを作成．oldevalなので切片だけで良い．
        resultList = []
        plt.subplot(subplot[0], subplot[1], index + 1)
        for slopeIndex, slope in enumerate(slopeList):
            for dosePoint, dose in enumerate(doseList[slopeIndex]):
                resultList.append([slope, (dosePoint + 1) * 10, checkEpsilon(name, dose)])
        data = pd.DataFrame(resultList, columns=["Slope", "MidPoint", "result"])
        heatmapData = pd.pivot_table(data=data, columns="Slope", index="MidPoint", values="result")
        ax = sns.heatmap(heatmapData, vmin=-1, vmax=1, cmap=cmap, linewidths=.3)
        ax.set_ylabel("MidPoint", fontsize=axizFontSize) # y軸
        ax.set_xlabel("Slope", fontsize=axizFontSize) # x軸
        ax.set_title("{} vs {}".format(name[0][:3], name[1][:3]), fontsize=titleFontSize)
        ax.tick_params(labelsize=labelSize)
        data.to_csv("{}/{}_{}.csv".format(csvdir, name[0], name[1]), index=False)
    fig.tight_layout()
    return fig

def oldeval_usecsv(dNameList, subplot, titleFontSize=16, axizFontSize=16, labelSize=16, csvdir="."):
    # 論文の評価関数を使用してヒートマップを作成する関数(作成したデータを使う場合)
    figsize = (subplot[1] * 10, subplot[0] * 10)
    if figsize: fig = plt.figure(figsize=figsize)
    else: fig = plt.figure()
    for index, name in enumerate(dNameList):
        plt.subplot(subplot[0], subplot[1], index + 1)
        data = pd.read_csv("{}/{}_{}.csv".format(csvdir, name[0], name[1]))
        heatmapData = pd.pivot_table(data=data, columns="Slope", index="MidPoint", values="result")
        ax = sns.heatmap(heatmapData, vmin=-1, vmax=1, cmap=cmap, linewidths=.3)
        ax.set_ylabel("MidPoint", fontsize=axizFontSize) # y軸
        ax.set_xlabel("Slope", fontsize=axizFontSize) # x軸
        ax.set_title("{} vs {}".format(name[0][:3], name[1][:3]), fontsize=titleFontSize)
        ax.tick_params(labelsize=labelSize)
    fig.tight_layout()
    return fig

def heatmap(dNameList, IC30, subplot, csvdir=".", axizFontSize=16, labelSize=16, splitNum=11):
    fig = plt.figure(figsize=(subplot[1] * 10, subplot[0] * 10))
    for index, name in enumerate(dNameList):
        plt.subplot(subplot[0], subplot[1], index + 1)
        drug = [makeDrugDatas(name[0]), makeDrugDatas(name[1])]
        doses = [[x, y] for x in np.linspace(0, IC30[name[0]] * 2, splitNum) for y in np.linspace(0, IC30[name[1]] * 2, splitNum)]
        resultList = []
        for dose in doses:
            resultList.append([round(dose[0], 2), round(dose[1], 2), sim(drug, dose)])
        data = pd.DataFrame(resultList, columns=["a1", "a2", "growth"])
        createHeatmap(data, name, axizFontSize=axizFontSize, labelSize=labelSize)
        data.to_csv("{}/{}_{}.csv".format(csvdir, name[0], name[1]), index=False)
    fig.tight_layout()

    return fig

def heatmap_usecsv(dNameList, subplot, csvdir=".", axizFontSize=16, labelSize=16):
    fig = plt.figure(figsize=(subplot[1] * 10, subplot[0] * 10))
    for index, name in enumerate(dNameList):
        plt.subplot(subplot[0], subplot[1], index + 1)
        data = pd.read_csv("{}/{}_{}.csv".format(csvdir, name[0], name[1]))
        createHeatmap(data, name, axizFontSize=axizFontSize, labelSize=labelSize)
    fig.tight_layout()

    return fig

if __name__ == "__main__":
    # 保存用ディレクトリの作成
    csvdir = "./results/ribo4/csv/virtual_drug"
    imgdir = "./results/ribo4/images/virtual_drug"
    makedir(csvdir)
    makedir(imgdir)

    ## drug data
    dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    # a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}

    # IC30 = calcIC(dNames, {name: 100 for name in dNames}, .3) # IC30を計算
    IC30 = {'Kanamycin': 0.6761401891708374, 'Streptmycin': 1.4652013778686523, 'Chloramphenicol': 21.09375, 'Tetracycline': 5.2734375}
    cmap = makeCmap() # 論文のデータと同じカラーマップ

    # dNameList = [[name, name] for name in dNames] # 同じ薬剤を２剤投与した場合．
    dNameList = itr.combinations(dNames, 2) # 異なる薬剤を２剤投与した場合．
    slopeList = [1./4, 1./2, 1., 2., 4.] # 傾きのリスト

    imgname = "{}/virtual_drug.png".format(imgdir)
    fig = heatmap(dNameList, IC30, (2, 3), csvdir, axizFontSize=40, labelSize=30)
    # fig = heatmap_usecsv(dNameList, (2, 3), csvdir, axizFontSize=40, labelSize=30)
    fig.savefig(imgname, dpi=300)
