# coding: utf-8
import matplotlib
matplotlib.use("Agg") # GUIのない環境でのおまじない
import sys
import pandas as pd
import itertools as itr
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np

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
                                     index = "a1",
                                     columns = "a2")
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
    a1DoseList = list(set(data["a1"].tolist()))[::20] # y軸に表示したい数のリスト
    a2DoseList = list(set(data["a2"].tolist()))[::20] # x軸に表示したい数のリスト

    # yticksの設定
    dataLenY = len(list(set(data["a1"].tolist()))) 
    ax.set_yticks(list(np.linspace(0.5, dataLenY - 0.5, len(a1DoseList))))
    ax.set_yticklabels(list(map(lambda x:"{:.3f}".format(x), a1DoseList)))
    
    # xticksの設定
    dataLenX = len(list(set(data["a2"].tolist()))) 
    ax.set_xticks(list(np.linspace(0.5, dataLenX - 0.5, len(a2DoseList))))
    ax.set_xticklabels(list(map(lambda x:"{:.3f}".format(x), a2DoseList)))

if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    csvdir = "results/ribo5/csv/sim100"
    
    # combinatorial simulation
    ## merge DataFiles
    # for modif in range(1, 3):
    #     for drugName in drugNameList:
    #         dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(drugName))
    #         fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
    #         df = mergeResults(fileNameList)
    #         df.to_csv("{}/modif{}/{}_merge.csv".format(csvdir, modif, "_".join(drugName)), index=False)

    ## SameDrug combination
    # drugNameList = [[name, name] for name in drugNames]
    # for modif in range(1, 3):
    #     dataLists = ["results/ribo5/csv/sim100/modif{}/{}_merge.csv".format(modif, "_".join(i)) for i in drugNameList]
    #     saveName = "results/ribo5/images/sameDrug_sim100_modif{}.png".format(modif)
    #     createHeatmap(drugNameList, dataLists, [2, 2], saveName)


    ## differentDrug combination
    # drugNameList = list(itr.combinations(drugNames, 2))
    # for modif in range(1, 3):
    #     dataLists = ["results/ribo5/csv/sim100/modif{}/{}_merge.csv".format(modif, "_".join(i)) for i in drugNameList]
    #     saveName = "results/ribo4/images/diffDrug_sim100_modif{}.png".format(modif)
    #     createHeatmap(drugNameList, dataLists, [3, 2], saveName)

    # combinatorial simulation(virtual drug)
    ## merge DataFiles
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["30s", "30s"], ["30s", "50s"]]
    nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for drugName in drugNameList for target in targetList]
    for modif in range(1, 3):
        csvdir = "results/ribo5/csv/sim100_v/modif{}".format(modif)
        for name in nameList:
            dirName = "{}/{}".format(csvdir, "_".join(name))
            fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
            df = mergeResults(fileNameList)
            df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(name)), index=False)
            
    ## create Image
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["30s", "30s"], ["30s", "50s"]]
    nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for target in targetList for drugName in drugNameList]
    for modif in range(1, 3):
        dataLists = ["results/ribo5/csv/sim100_v/modif{}/{}_merge.csv".format(modif, "_".join(name)) for name in nameList]
        saveName = "results/ribo4/images/virtualDrug_sim100_modif{}.png".format(modif)
        createHeatmap(nameList, dataLists, [3, 2], saveName)



