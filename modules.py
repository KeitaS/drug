# coding: utf-8

from ribo7 import *
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools
import math


def generate_cmap(colors):
    """自分で定義したカラーマップを返す(線形補完)"""
    from matplotlib.colors import LinearSegmentedColormap
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)

def makedir(dirname):
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    del(os)

def sim(drugs, dose):
    for i in range(len(drugs)):
        drugs[i]["dose"] = dose[i]
    result, legend = run(drugs, legend=["r_u"])
    return calcGrowthRate(result[-1][1])

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

def checkEpsilon(dName, dose): # nature genesis 2006's evaluation
    x = sim([makeDrugDatas(dName[0])], [dose[0]])
    y = sim([makeDrugDatas(dName[1])], [dose[1]])
    xy = sim([makeDrugDatas(dName[0]), makeDrugDatas(dName[1])], dose)
    return (xy - x * y) / abs(min(x, y) - x * y)

def growthPlot(dNames, IC30, savename=""): # glowth calc.
    for index, dName in enumerate(dNames):
        drugs = [makeDrugDatas(dName)]
        doses = np.linspace(0, IC30[dName] * 2, 101) # 0 - IC30 * 2
        x = math.ceil(math.sqrt(len(dNames)))
        plt.subplot(int(math.ceil(len(dNames) / x)), int(x), index + 1) # subplotの動的確保
        plt.plot(doses, [sim(drugs, [dose]) for dose in doses])
        plt.title(dName)
        plt.xlabel("Dose")
        plt.ylabel("Growth Rate ($\lambda / \lambda_0$)")
    plt.tight_layout()
    if savename: plt.savefig(savename, dpi=300)
    else: plt.show()


def divideFigure(drugNames):
    """
        drugNamesに応じて，figureの分割数を計算するモジュール．
    """
    x = int(math.ceil(math.sqrt(len(drugNames))))
    y = int(math.ceil(len(drugNames) / x))
    return (x, y)


def createHeatmap(data, drugNames, cbar=False, cmap=False):
    """
        function of create Heatmap.
        data: pandas data.
        cmap: color map of heatmap.
    """
    if not cmap: cmap = sns.diverging_palette(220, 10, as_cmap=True) # coler map
    heatmap = pd.pivot_table(data=data, values="growth", index="a1", columns="a2") # heatmap data
    sns.heatmap(heatmap, cbar=False, cmap=cmap) # create Heatmap
    plt.ylabel(drugNames[0], fontsize=16) # create ylabel
    plt.xlabel(drugNames[1], fontsize=16) # create xlabel

    ax = plt.gca()
    ax.invert_yaxis() # reverse ylabel


def heatmap(drugNames, IC30, split=11, savename="", csvdir=""):
    x, y = divideFigure(drugNames) # calc division number
    plt.figure(figsize=(x * 5, y * 5)) # set figure size

    for index, drugNameList in enumerate(drugNames):
        drugs = [makeDrugDatas(drugNameList[0]), makeDrugDatas(drugNameList[1])] # create drug data
        doses = [[x, y] for x in np.linspace(0, IC30[drugNameList[0]] * 2, split) for y in np.linspace(0, IC30[drugNameList[1]] * 2, split)] # create dose patern
        data = pd.DataFrame([[round(dose[0], 2), round(dose[1], 2), sim(drugs, dose)] for dose in doses],
                            columns=["a1", "a2", "growth"]) # simulation & create data

        if csvdir: data.to_csv("{}/{}_{}.csv".format(csvdir, drugNameList[0], drugNameList[1]), index=False) # make csv data

        plt.subplot(y, x, index + 1) # create subplot area
        createHeatmap(data, drugNameList)

    plt.tight_layout()
    if savename: plt.savefig(savename, dpi=300)
    else: plt.show()

def neweval(dNames, IC30, split=11, savename="", figsize=(12, 9), csvdir=""):
    plt.figure(figsize=figsize)
    cmap = generate_cmap(["mediumblue", "white", "orangered"])
    slopeList = [1./4, 1./2, 1., 2., 4.] # 傾き
    divnum = 11 # lineを見るときのdivision number

    for index, dNameList in enumerate(dNames):
        drugs = [makeDrugDatas(dNameList[0]), makeDrugDatas(dNameList[1])]
        # 中点をリスト化
        midPointList = [x for x in itertools.zip_longest(np.linspace(0, IC30[dNameList[0]], split), np.linspace(0, IC30[dNameList[1]], split))][1:]
        # slopeに応じたdoselistの生成
        dosesList = [[[x for x in itertools.zip_longest(np.linspace(0, midPoint[0] * (1 + slope), 11), np.linspace(0, midPoint[1] * (1 + (1 / slope)), 11)[::-1])] for midPoint in midPointList] for slope in slopeList]
        resultList = []
        for slope, a in enumerate(dosesList):
            for midPoint, doses in enumerate(a):
                resultList.append([slopeList[slope], midPoint + 1, [sim(drugs, dose) for dose in doses]])

        data = pd.DataFrame(resultList, columns=["Slope", "MidPoint", "GrowthList"])
        data["LinerType"] = [checkLinerType(inp) for inp in data["GrowthList"]]
        if csvdir: data.to_csv("{}/{}_{}.csv".format(csvdir, dNameList[0][:3], dNameList[1][:3]), index=False)
        heatmap = pd.pivot_table(data=data, values="LinerType", index="MidPoint", columns="Slope")
        # make subplot
        x = math.ceil(math.sqrt(len(dNames)))
        plt.subplot(int(math.ceil(len(dNames) / x)), int(x), index + 1)
        # heatmap
        sns.heatmap(heatmap, vmin=-1, vmax=1, cmap=cmap, linewidths=.3, annot=True, annot_kws={"size": 7})
        plt.gca().invert_yaxis()
        plt.ylabel("Midpoint")
        plt.xlabel("Slope")
        plt.title("{} vs {}".format(dNameList[0][:3], dNameList[1][:3]))
    plt.tight_layout()
    plt.tick_params(labelsize=10)
    if savename: plt.savefig(savename, dpi=300)
    else: plt.show()

def oldeval(dNames, IC30, split=11, savename="", figsize=(12, 9)):
    plt.figure(figsize=figsize)
    cmap = generate_cmap(["mediumblue", "white", "orangered"])
    slopeList = [1./4, 1./2, 1., 2., 4.] # 傾き
    divnum = 11 # lineを見るときのdivision number

    for index, dNameList in enumerate(dNames):
        drugs = [makeDrugDatas(dNameList[0]), makeDrugDatas(dNameList[1])]
        midPointList = [x for x in itertools.zip_longest(np.linspace(0, IC30[dNameList[0]] * 2, split), np.linspace(0, IC30[dNameList[1]] * 2, split))][1:]
        dosesList = [[[x for x in itertools.zip_longest(np.linspace(0, midPoint[0] * (1 + slope), 11), np.linspace(0, midPoint[1] * (1 + (1 / slope)), 11)[::-1])] for midPoint in midPointList] for slope in slopeList]
        resultList = []
        for slope, a in enumerate(dosesList):
            for midPoint, doses in enumerate(a):
                resultList.append([slopeList[slope], midPoint + 1, checkEpsilon([sim(drugs, dose) for dose in doses])])
        data = pd.DataFrame(resultList, columns=["Slope", "Midpoint", "result"])
        heatmap = pd.pivot_table(data=data, values="result", index="Midpoint", columns="Slope")
        x = math.ceil(math.sqrt(len(dNames)))
        plt.subplot(int(math.ceil(len(dNames) / x)), int(x), index + 1)
        sns.heatmap(heatmap, vmin=-1, vmax=1, cmap=cmap, linewidths=.3, annot=True, annot_kws={"size": 7})
        plt.ylabel("Midpoint")
        plt.xlabel("Slope")
        plt.title("{} vs {}".format(dNameList[0][:3], dNameList[1][:3]))
    plt.tight_layout()
    plt.tick_params(labelsize=10)
    if savename: plt.savefig(savename, dpi=300)
    else: plt.show()

def diffTargetCheck(dNames, IC30, figsize=(12, 9), savename=""):
    """
        同薬剤のターゲットを変えて組み合わせ投与する．
    """
    plt.figure(figsize=figsize)
    for index, dName in enumerate(dNames):
        drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
        if drugs[0]["target"] == "r30":
            drugs[1]["target"] = "r50"
        else:
            drugs[1]["target"] = "r30"
        doses = [[x, y] for x in np.linspace(0, IC30[dName] * 2, 11) for y in np.linspace(0, IC30[dName] * 2, 11)]
        data = pd.DataFrame([[round(dose[0], 2), round(dose[1], 2), sim(drugs, dose)] for dose in doses], columns=[drugs[0]["target"], drugs[1]["target"], "growth"])

        heatmap = pd.pivot_table(data=data, values="growth", index=drugs[0]["target"], columns=drugs[1]["target"])
        x = math.ceil(math.sqrt(len(dNames)))
        plt.subplot(int(math.ceil(len(dNames) / x)), int(x), index + 1)
        sns.heatmap(heatmap)
        plt.ylabel(drugs[0]["target"])
        plt.xlabel(drugs[1]["target"])
        plt.title(dName)
    plt.tight_layout()
    plt.savefig(savename, dpi=300)


# テスト用モジュール
def sim_ribo5(drugs, dose, modif=0):
    """
        ribo5を使ってsim．
    """
    import ribo5
    for i in range(len(drugs)):
        drugs[i]["dose"] = dose[i]
    result, legend = ribo5.run(drugs, legend=["r_u"], inpData={"modif": modif})
    return calcGrowthRate(result[-1][1])


def ribo5_targetChange(drugNames, IC30, target=[], modif=0):
    import ribo5
    drugs = [ribo5.makeDrugDatas(drugNames[0]), ribo5.makeDrugDatas(drugNames[1])]
    if target:
        drugs[0]["type"] = target[0]
        drugs[1]["type"] = target[1]
    doses = [[x, y] for x in np.linspace(0, IC30[drugNames[0]] * 2, 11) for y in np.linspace(0, IC30[drugNames[1]] * 2, 11)]
    data = pd.DataFrame([[round(dose[0], 2), round(dose[1], 2), sim_ribo5(drugs, dose, modif=modif)] for dose in doses], columns=["a1", "a2", "growth"])
    return data

def ribo5_heatmap(drugNames, IC30, split=11, savename="", csvdir=""):
    import ribo5
    x, y = divideFigure(drugNames) # calc division number
    plt.figure(figsize=(x * 5, y * 5)) # set figure size

    for index, drugNameList in enumerate(drugNames):
        drugs = [ribo5.makeDrugDatas(drugNameList[0]), ribo5.makeDrugDatas(drugNameList[1])] # create drug data
        doses = [[x, y] for x in np.linspace(0, IC30[drugNameList[0]] * 2, split) for y in np.linspace(0, IC30[drugNameList[1]] * 2, split)] # create dose patern
        data = pd.DataFrame([[round(dose[0], 2), round(dose[1], 2), sim_ribo5(drugs, dose)] for dose in doses],
                            columns=["a1", "a2", "growth"]) # simulation & create data

        if csvdir: data.to_csv("{}/{}_{}.csv".format(csvdir, drugNameList[0], drugNameList[1]), index=False) # make csv data

        plt.subplot(y, x, index + 1) # create subplot area
        createHeatmap(data, drugNameList)

    plt.tight_layout()
    if savename: plt.savefig(savename, dpi=300)
    else: plt.show()
