# coding: utf-8
import matplotlib
matplotlib.use("Agg") # GUIのない環境でのおまじない
import sys
import pandas as pd
import itertools as itr
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
from ribo4 import makeCmap
from ribo4 import generate_cmap

def mergeResults(fileNameList):
    """
        並列計算した結果を１つのファイルにまとめる関数
        fileNameList : ファイル名をリストにしたもの．
    """
    df = pd.DataFrame()
    for fileName in fileNameList:
        data = pd.read_csv(fileName)
        df = pd.concat([df, data])
    df = df.reset_index(drop=True)
    return df

def createHeatmap(drugNames, dataLists, subplot, saveName, xy=["a2", "a1"], simType="growth"):
    """
        drugNames : csvのファイル名で使用している薬剤名のリスト
        csvdir    : csvが保存されているディレクトリ名
        subplot   : グラフの数を行列で[横, 縦]
        saveName  : 保存するときの名前
        xy        : 軸にするDataFrameのcolumnsの名前（x軸，y軸の順)
        simType   : シミュレーションのタイプ（growth, epsilon, LynerType)
    """
    fig = plt.figure(figsize=(subplot[0] * 100 / 9, subplot[1] * 10))
    cbar_ax = fig.add_axes([.92, .1, .02, .8])
    cmapStatus = {"growth" : {"cmap": sns.diverging_palette(220, 10, as_cmap=True), "vmax": 1., "vmin": 0.},
                  "growthRate" : {"cmap": sns.diverging_palette(220, 10, as_cmap=True), "vmax": 1., "vmin": 0.},
                  "epsilon": {"cmap": makeCmap(), "vmax": 1., "vmin": -1.}, 
                  "LinerType": {"cmap": generate_cmap(), "vmax": .3, "vmin": -.3}}
    for index, drugName in enumerate(drugNames):
        plt.subplot(subplot[1], subplot[0], index + 1)
        data = pd.read_csv(dataLists[index])
        heatmapData = pd.pivot_table(data = data,
                                     values = simType,
                                     index = xy[1],
                                     columns = xy[0])
        ax = sns.heatmap(heatmapData,
                         annot = simType == "LinerType",
                         annot_kws={"size": 20},
                         fmt="1.3f",
                         cbar = index == 0,
                         cmap = cmapStatus[simType]["cmap"],
                         vmax = cmapStatus[simType]["vmax"],
                         vmin = cmapStatus[simType]["vmin"],
                         cbar_ax = None if index else cbar_ax,
                         square = simType != "LinerType")

        ax.invert_yaxis() # y軸の上下を変える
        if simType != "LinerType":
            ax.set_ylabel(drugName[0], fontsize=30) # y軸のラベル
            ax.set_xlabel(drugName[1], fontsize=30) # x軸のラベル
            setTickLabel(data, ax, xy) # 軸の設定
        else :
            ax.set_ylabel(xy[1], fontsize=30) 
            ax.set_xlabel(xy[0], fontsize=30) # x軸のラベル
            ax.set_title("{} vs {}".format(drugName[0], drugName[1]), fontsize=30)

        ax.tick_params(labelsize=24)

    cbar_ax.tick_params(labelsize = 30)
    fig.tight_layout(rect = [0, 0, .9, 1.])
    plt.savefig(saveName, dpi=300)


def setTickLabel(data, ax, xy):
    a1DoseList = sorted(list(set(data[xy[1]].tolist())))[::20] # y軸に表示したい数のリスト
    a2DoseList = sorted(list(set(data[xy[0]].tolist())))[::20] # x軸に表示したい数のリスト

    # yticksの設定
    dataLenY = len(list(set(data[xy[1]].tolist()))) 
    ax.set_yticks(list(np.linspace(0.5, dataLenY - 0.5, len(a1DoseList))))
    ax.set_yticklabels(list(map(lambda x:"{:.3f}".format(x), a1DoseList)))
    
    # xticksの設定
    dataLenX = len(list(set(data[xy[0]].tolist()))) 
    ax.set_xticks(list(np.linspace(0.5, dataLenX - 0.5, len(a2DoseList))))
    ax.set_xticklabels(list(map(lambda x:"{:.3f}".format(x), a2DoseList)))

def makeCmap(c_range={"red": 30, "pink": 5, "white": 10, "light_green": 5, "green": 13, "blue": 17},
            c_num = {"red": "#ff0000", "pink": "#ffc0cb", "white": "#ffffff", "light_green": "#90ee90", "green": "#008000", "blue": "#0000ff"},
            c_list = ["red", "pink", "white", "light_green", "green", "blue"]): # eをかませるためのカラーマップを作る関数
    """自分で定義したカラーマップを返す(固定)"""

    c_result = []
    for color in c_list:
        for i in range(c_range[color]):
            c_result.append(c_num[color])
    cmap = matplotlib.colors.ListedColormap(c_result)
    return cmap

def generate_cmap(colors=["mediumblue", "white", "orangered"]):
    """自分で定義したカラーマップを返す(線形補完)"""
    from matplotlib.colors import LinearSegmentedColormap
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


