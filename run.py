# coding: utf-8
from ribo7 import *
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools
import math

dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
# IC30 = calcIC(dNames, {x: 50 for x in dNames}, 0.3)
IC30 = {'Kanamycin': 0.4598379135131836, 'Tetracycline': 0.732421875, 'Chloramphenicol': 4.3212890625, 'Streptmycin': 0.9600967168807983}

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

def checkEpsilon(dName, dose):
    x = sim([makeDrugDatas(dName[0])], [dose[0]])
    y = sim([makeDrugDatas(dName[1])], [dose[1]])
    xy = sim([makeDrugDatas(dName[0]), makeDrugDatas(dName[1])], dose)
    return (xy - x * y) / abs(min(x, y) - x * y)

def growthPlot(dNames, IC30, type=0, savename=""):
    for index, dName in enumerate(dNames):
        drugs = [makeDrugDatas(dName)]
        doses = np.linspace(0, IC30[dName] * 2, 101)
        plt.subplot(2, 2, index + 1)
        plt.plot(doses, [sim(drugs, [dose]) for dose in doses])
        plt.title(dName)
        plt.xlabel("Dose")
        plt.ylabel("Growth Rate ($\lambda / \lambda_0$)")
    plt.tight_layout()
    if type == 0: plt.show()
    elif type == 1: plt.savefig(savename, dpi=300)

def heatmap(dNames, IC30, split=11, type=0, savename="", figsize=(12, 9)):
    plt.figure(figsize=figsize)
    for index, dNameList in enumerate(dNames):
        drugs = [makeDrugDatas(dNameList[0]), makeDrugDatas(dNameList[1])]
        doses = [[x, y] for x in np.linspace(0, IC30[dNameList[0]] * 2, split) for y in np.linspace(0, IC30[dNameList[1]] * 2, split)]
        data = pd.DataFrame([[round(dose[0], 2), round(dose[1], 2), sim(drugs, dose)] for dose in doses], columns=["a1", "a2", "growth"])
        heatmap = pd.pivot_table(data=data, values="growth", index="a1", columns="a2")
        x = math.ceil(math.sqrt(len(dNames)))
        plt.subplot(int(math.ceil(len(dNames) / x)), int(x), index + 1)
        sns.heatmap(heatmap)
        plt.ylabel(dNameList[0])
        plt.xlabel(dNameList[1])
    plt.tight_layout()
    plt.tick_params(labelsize=10)
    if type == 0: plt.show()
    elif type == 1: plt.savefig(savename, dpi=300)

def neweval(dNames, IC30, split=11, type=0, savename="", figsize=(12, 9)):
    plt.figure(figsize=figsize)
    cmap = generate_cmap(["mediumblue", "white", "orangered"])
    slopeList = [1./4, 1./2, 1., 2., 4.] # 傾き
    divnum = 11 # lineを見るときのdivision number

    for index, dNameList in enumerate(dNames):
        drugs = [makeDrugDatas(dNameList[0]), makeDrugDatas(dNameList[1])]
        midPointList = [x for x in itertools.zip_longest(np.linspace(0, IC30[dNameList[0]], split), np.linspace(0, IC30[dNameList[1]], split))][1:]
        dosesList = [[[x for x in itertools.zip_longest(np.linspace(0, midPoint[0] * (1 + slope), 11), np.linspace(0, midPoint[1] * (1 + (1 / slope)), 11)[::-1])] for midPoint in midPointList] for slope in slopeList]
        resultList = []
        for slope, a in enumerate(dosesList):
            for midPoint, doses in enumerate(a):
                resultList.append([slopeList[slope], midPoint + 1, checkLinerType([sim(drugs, dose) for dose in doses])])
        data = pd.DataFrame(resultList, columns=["Slope", "Midpoint", "result"])

        heatmap = pd.pivot_table(data=data, values="result", index="Midpoint", columns="Slope")
        x = math.ceil(math.sqrt(len(dNames)))
        plt.subplot(int(math.ceil(len(dNames) / x)), int(x), index + 1)
        sns.heatmap(heatmap, vmin=-1, vmax=1, cmap=cmap, linewidths=.3, annot=True, annot_kws={"size": 7})
        plt.gca().invert_yaxis()
        plt.ylabel("Midpoint")
        plt.xlabel("Slope")
        plt.title("{} vs {}".format(dNameList[0][:3], dNameList[1][:3]))
    plt.tight_layout()
    plt.tick_params(labelsize=10)
    if type == 0: plt.show()
    elif type == 1: plt.savefig(savename, dpi=300)

def oldeval(dNames, IC30, split=11, type=0, savename="", figsize=(12, 9)):
    plt.figure(figsize=figsize)
    cmap = generate_cmap(["mediumblue", "white", "orangered"])
    slopeList = [1./4, 1./2, 1., 2., 4.] # 傾き
    divnum = 11 # lineを見るときのdivision number

    for index, dNameList in enumerate(dNames):
        drugs = [makeDrugDatas(dNameList[0]), makeDrugDatas(dNameList[1])]
        midPointList = [x for x in itertools.zip_longest(np.linspace(0, IC30[dNameList[0]], split), np.linspace(0, IC30[dNameList[1]], split))][1:]
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
    if type == 0: plt.show()
    elif type == 1: plt.savefig(savename, dpi=300)

#
# 保存するディレクトリの設定
savedir = "images/ribo7"
makedir(savedir)
# growthPlot(dNames, IC30, 1, savedir+"/single_all_growthrate.png")
# heatmap(list(itertools.combinations(dNames, 2)), IC30, 11, 1, savedir + "/double_heatmap_comb.png", figsize=(18, 9))
# heatmap([[dName, dName] for dName in dNames], IC30, 11, 1, savedir + "/double_heatmap.png", figsize=(12, 9))
neweval(list(itertools.combinations(dNames, 2)), IC30, 11, 1, savedir + "/double_neweval_comb.png", figsize=(18, 9))
neweval([[dName, dName] for dName in dNames], IC30, 11, 1, savedir + "/double_neweval.png", figsize=(12, 9))
