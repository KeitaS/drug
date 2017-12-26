#coding: utf-8
import matplotlib
matplotlib.use("Agg") # guiのない環境でのおまじない
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools as itr
import numpy as np

def makedir(dirname):
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    del(os)

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

def createSingleGraph(dataPath, savename):
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    plt.figure(figsize=(16, 12))
    for index, drugName in enumerate(drugNames):
        data = pd.read_csv("{}/{}.csv".format(dataPath, drugName))
        plt.subplot(2, 2, index + 1)
        plt.plot(data["dose"], data["growth"])
        plt.xlabel("Dose", fontsize=20)
        plt.ylabel("Growth Rate", fontsize=20)
        plt.title(drugName, fontsize=30)
        plt.tick_params(labelsize=14)
    plt.tight_layout()
    plt.savefig(savename, dpi=300)

def createHeatmap(drugNames, dataLists, subplot, saveName):
    """
        drugNames : csvのファイル名で使用している薬剤名のリスト
        csvdir    : csvが保存されているディレクトリ名
        subplot   : グラフの数を行列で[横, 縦]
    """
    fig = plt.figure(figsize=(subplot[0] * 100 / 9, subplot[1] * 10))
    cbar_ax = fig.add_axes([.92, .1, .02, .8])
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    for index, drugName in enumerate(drugNames):
        plt.subplot(subplot[1], subplot[0], index + 1)
        data = pd.read_csv(dataLists[index])
        heatmapData = pd.pivot_table(data = data,
                                     values = "growth",
                                     index = "dose1",
                                     columns = "dose2")
        ax = sns.heatmap(heatmapData,
                         cbar = index == 0,
                         cmap = cmap,
                         cbar_ax = None if index else cbar_ax,
                         square = True)

        ax.invert_yaxis() # y軸の上下を変える
        setTickLabel(data, ax) # 軸の設定
        ax.set_ylabel(drugName[0], fontsize = 30) # y軸のラベル
        ax.set_xlabel(drugName[1], fontsize = 30) # x軸のラベル
        ax.tick_params(labelsize=24)

    cbar_ax.tick_params(labelsize = 30)
    fig.tight_layout(rect = [0, 0, .9, 1.])
    plt.savefig(saveName, dpi=300)


def setTickLabel(data, ax):
    dose1DoseList = list(set(data["dose1"].tolist()))[::20] # y軸に表示したい数のリスト
    dose2DoseList = list(set(data["dose2"].tolist()))[::20] # x軸に表示したい数のリスト

    # yticksの設定
    dataLenY = len(list(set(data["dose1"].tolist()))) 
    ax.set_yticks(list(np.linspace(0.5, dataLenY - 0.5, len(dose1DoseList))))
    ax.set_yticklabels(list(map(lambda x:"{:.3f}".format(x), dose1DoseList)))
    
    # xticksの設定
    dataLenX = len(list(set(data["dose2"].tolist()))) 
    ax.set_xticks(list(np.linspace(0.5, dataLenX - 0.5, len(dose2DoseList))))
    ax.set_xticklabels(list(map(lambda x:"{:.3f}".format(x), dose2DoseList)))


if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    

    
    # single simulation 
    singleDataPath = "results/ribo8/single/csv"
    imagesDirPath = "results/ribo8/images"
    makedir(singleDataPath)
    makedir(imagesDirPath)

    # combinatorial simulation
    ## merge DataFiles
    # drugNameList = itr.combinations_with_replacement(drugNames, 2)
    # csvdir = "results/ribo8/double/normal"

    # for drugName in drugNameList:
    #     dirName = "{}/{}".format(csvdir, "_".join(drugName))
    #     fileNameList = ["{}/{}_{}.csv".format(dirName, "_".join(drugName), num) for num in range(101)]
    #     df = mergeResults(fileNameList)
    #     df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(drugName)), index=False)

    ## SameDrug combination
    drugNameList = [[name, name] for name in drugNames]
    dataLists = ["results/ribo8/double/normal/{}_merge.csv".format("_".join(i)) for i in drugNameList]
    doubleSaveName = "results/ribo8/images/sameDrug.png"
    createHeatmap(drugNameList, dataLists, [2, 2], doubleSaveName)

    ## differentDrug combination
    drugNameList = list(itr.combinations(drugNames, 2))
    dataLists = ["results/ribo8/double/normal/{}_merge.csv".format("_".join(i)) for i in drugNameList]
    doubleSaveName = "results/ribo8/images/diffDrug.png"
    createHeatmap(drugNameList, dataLists, [3, 2], doubleSaveName)
