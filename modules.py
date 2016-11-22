# coding: utf-8
from ribo5 import *
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pylab as plt


def makedir(dirname):
    """ディレクトリ作成関数"""
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    del(os)

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


def epsilon(x, y, val):
    """
    Nature genetics 2006のepsilonを計算する関数
    e = (Wxy - WxWy)/|min(Wx, Wy) - WxWy|
    """
    result = (val - x * y) / abs(min(x, y) - x * y)
    return result


def doseResponse(drugs, dose, inpData={"K_ma": 15., "modif": 1}):
    """
        run して，Growth rateを返す関数
    """
    for i in range(len(drugs)):
        drugs[i]["dose"] = dose[i]
    result, legend = run(drugs, step=100, inpData=inpData, legend=["r_u"])
    result = calcGrowthrate(result[-1][1])
    return result

def createSlopedose(slope, midPointList, divnum=11):
    """
        midPointに応じて，該当するSlopeで作成したDoseの組み合わせリストを返す関数
    """
    doseX = np.linspace(0, midPointList[0] * (1 + slope), divnum)
    doseY = np.linspace(0, midPointList[1] * (1 + (1 / slope)), divnum)[::-1]
    return_list = [[doseX[i], doseY[i]] for i in range(len(doseX))]
    return [[doseX[i], doseY[i]] for i in range(len(doseX))]

def calcEpsilon(dNames, doses, inpData={"K_ma": 15, "modif": 1}):
    """
    """
    result_list = []
    for i in range(3):
        if i < 2:
            drugs = [makeDrugDatas(dNames[i])]
            result_list.append(doseResponse(drugs, [doses[i]], inpData=inpData))
        else:
            drugs = [makeDrugDatas(dNames[0]), makeDrugDatas(dNames[1])]
            result_list.append(doseResponse(drugs, doses, inpData=inpData))

    return epsilon(result_list[0], result_list[1], result_list[2])




def growthHeatmap(data, values, index, columns, title="", xlabel="", ylabel=""):
    heatmap = pd.pivot_table(data=data, values=values , index=index, columns=columns)
    sns.heatmap(heatmap)

    # ラベルの設定
    if ylabel: plt.ylabel(ylabel)
    else: plt.ylabel(index)

    if xlabel: plt.xlabel(xlabel)
    else: plt.xlabel(columns)

    # タイトルの設定
    if title:
        plt.title(title)

    # ラベルの文字サイズ
    plt.tick_params(labelsize=7)

def evalHeatmap(data, cmap, values, index, columns, title="", xlabel="", ylabel=""):
    heatmap = pd.pivot_table(data=data, values=values , index=index, columns=columns)
    sns.heatmap(heatmap, vmin=-1, vmax=1, cmap=cmap, linewidths=.3, annot=True, annot_kws={"size": 7})


    # ラベルの設定
    if ylabel: plt.ylabel(ylabel)
    else: plt.ylabel(index)

    if xlabel: plt.xlabel(xlabel)
    else: plt.xlabel(columns)

    # タイトルの設定
    if title:
        plt.title(title)

    # ラベルの文字サイズ
    plt.tick_params(labelsize=7)
