# coding: utf-8
from ribo7 import *
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools as itr
import math
from modules import *

def neweval(dNameList, IC30, subplot, csvdir=".", titleFontSize=16, axizFontSize=16, labelSize=16, splitNum=11, cmap=generate_cmap(), flag=True):
    # 論文の評価関数を使用してヒートマップを作成する関数
    figsize = (subplot[1] * 10, subplot[0] * 10)
    fig = plt.figure(figsize=figsize)
    # fig, axn = plt.subplots(subplot[0], subplot[1])
    # cbar_ax = fig.add_axes([.91, .3, .03, .4])
    slopeList = [1./4., 1./2., 1., 2., 4.]
    for index, name in enumerate(dNameList):
        drugs = [makeDrugDatas(name[0]), makeDrugDatas(name[1])]

        if flag:
            drugs[0]["target"] = "r30"
            drugs[1]["target"] = "r30"

        midPointList = [x for x in itr.zip_longest(np.linspace(0, IC30[name[0]], 11), np.linspace(0, IC30[name[1]], 11))][1:] # 0, 0を除いた10点．中点なので，IC30でOK．
        doseList = [[[x for x in itr.zip_longest(np.linspace(0, midPoint[0] * (1 + slope), 11), np.linspace(0, midPoint[1] * (1 + (1 / slope)), 11)[::-1])] for midPoint in midPointList] for slope in slopeList]
        resultList = []

        plt.subplot(subplot[0], subplot[1], index + 1)
        for slope, midPointList in enumerate(doseList):
            for midPoint, doses in enumerate(midPointList):
                resultList.append([(midPoint + 1) * 10, slopeList[slope], [sim(drugs, dose) for dose in doses]])
        data = pd.DataFrame(resultList, columns=["MidPoint", "Slope", "resultList"])
        data["LinerType"] = [checkLinerType(inp) for inp in data["resultList"]]
        heatmapData = pd.pivot_table(data=data, index="MidPoint", columns="Slope", values="LinerType") # valueをepsilonに
        # ax = sns.heatmap(heatmapData, vmin=-1, vmax=1, cmap=cmap,
        #                  cbar=index == 0, cbar_ax = None if index else cbar_ax,
        #                  linewidths=.2, annot=True, annot_kws={"size": 4}, fmt="1.2f")
        ax = sns.heatmap(heatmapData, vmin=-1, vmax=1, cmap=cmap,
                        cbar=False, linewidths=.2, annot=True, fmt="1.2f", annot_kws={"size": 30})
        ax.set_ylabel("MidPoint", fontsize=axizFontSize) # y軸
        ax.set_xlabel("Slope", fontsize=axizFontSize) # x軸
        ax.set_title("{} vs {}".format(name[0][:3], name[1][:3]), fontsize=titleFontSize)
        ax.tick_params(labelsize=labelSize)
        ax.invert_yaxis()

        data.to_csv("{}/{}_{}.csv".format(csvdir, name[0], name[1]), index=False)
    fig.tight_layout()
    return fig

def neweval_usecsv(dNameList, subplot, csvdir=".", titleFontSize=16, axizFontSize=16, labelSize=16, splitNum=11, cmap=generate_cmap()):
    # 論文の評価関数を使用してヒートマップを作成する関数
    fig = plt.figure(figsize=(subplot[1] * 10, subplot[0] * 10))
    # fig, axn = plt.subplots(subplot[0], subplot[1])
    # cbar_ax = fig.add_axes([.91, .35, .01, .3])
    for index, name in enumerate(dNameList):
        plt.subplot(subplot[0], subplot[1], index + 1)
        data = pd.read_csv("{}/{}_{}.csv".format(csvdir, name[0], name[1]))
        heatmapData = pd.pivot_table(data=data, index="MidPoint", columns="Slope", values="LinerType") # valueをepsilonに
        # ax = sns.heatmap(heatmapData, vmin=-1, vmax=1, cmap=cmap,
        #                  cbar=index == 0, cbar_ax = None if index else cbar_ax,
        #                  linewidths=.2, annot=True, annot_kws={"size": 4}, fmt="1.2f")
        ax = sns.heatmap(heatmapData, vmin=-1, vmax=1, cmap=cmap,
                        cbar=False, linewidths=.2, annot=True, fmt="1.2f", annot_kws={"size": 30})
        ax.set_ylabel("MidPoint", fontsize=axizFontSize) # y軸
        ax.set_xlabel("Slope", fontsize=axizFontSize) # x軸
        ax.set_title("{} vs {}".format(name[0][:3], name[1][:3]), fontsize=titleFontSize)
        ax.tick_params(labelsize=labelSize)
        ax.invert_yaxis()

    # fig.tight_layout(rect=[0, 0, .9 , 1])
    fig.tight_layout()
    return fig

def createHeatmap(data, drugNames, cbar=False, cmap=False, axizFontSize=16, labelSize=16):
    """
        function of create Heatmap.
        data: pandas data.
        cmap: color map of heatmap.
    """
    if not cmap: cmap = sns.diverging_palette(220, 10, as_cmap=True) # coler map
    heatmapData = pd.pivot_table(data=data, values="growth", index="a1", columns="a2") # heatmap data
    ax = sns.heatmap(heatmapData, cbar=cbar, cmap=cmap, square=True) # create Heatmap
    ax.invert_yaxis()

    setTickLabel(data, ax)

    ax.set_ylabel(drugNames[0], fontsize=axizFontSize) # create ylabel
    ax.set_xlabel(drugNames[1], fontsize=axizFontSize) # create xlabel

    ## virtual drug
    # if drugNames[0] == "Streptmycin": ax.set_ylabel("Pattern A", fontsize=axizFontSize) # create ylabel
    # else : ax.set_ylabel("Pattern B", fontsize=axizFontSize) # create ylabe
    # if drugNames    [1] == "Streptmycin": ax.set_xlabel("Pattern A", fontsize=axizFontSize) # create xlabel
    # else: ax.set_xlabel("Pattern B", fontsize=axizFontSize) # create xlabel

    ax.tick_params(labelsize=labelSize)

def heatmap(dNameList, IC30, subplot, csvdir=".", axizFontSize=16, labelSize=16, splitNum=11, flag=False):
    fig = plt.figure(figsize=(subplot[1] * 10, subplot[0] * 10))
    for index, name in enumerate(dNameList):
        plt.subplot(subplot[0], subplot[1], index + 1)
        drug = [makeDrugDatas(name[0]), makeDrugDatas(name[1])]

        if flag:
            drug[0]["target"] = "r30"
            drug[1]["target"] = "r30"

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

def setTickLabel(data, ax):
    a1DoseList = list(set(data["a1"].tolist()))[::2] # y軸に表示したい数のリスト
    a2DoseList = list(set(data["a2"].tolist()))[::2] # x軸に表示したい数のリスト

    ax.set_xticks(list(np.linspace(0.5, 10.5, len(a2DoseList)))) # xticksのせてい
    ax.set_xticklabels(list(map(str, a2DoseList)))
    ax.set_yticks(list(np.linspace(0.5, 10.5, len(a1DoseList))))
    ax.set_yticklabels(list(map(str, a1DoseList)))


dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
# IC30 = calcIC(dNames, {x: 50 for x in dNames}, 0.3)
# print(IC30)
IC30 = {'Streptmycin': 1.46942138671875, 'Kanamycin': 0.6775975227355957, 'Tetracycline': 5.2734375, 'Chloramphenicol': 21.09375}

csvdir = "results/ribo7/csv/neweval"
imgdir = "results/ribo7/images/neweval"

# dNameList = [[name, name] for name in dNames] # 同じ薬剤を２剤投与した場合．
dNameList = itr.combinations(dNames, 2) # 異なる薬剤を２剤投与した場合．
# dNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]

fig = neweval_usecsv(dNameList, (2,3), csvdir, axizFontSize=40, labelSize=30, titleFontSize=30)
fig.savefig("{}/combination_drug.png".format(imgdir), dpi=300)
