# coding: utf-8
import matplotlib
matplotlib.use("Agg") # GUI$B$N$J$$4D6-$G$N$*$^$8$J$$(B
import sys
import pandas as pd
import itertools as itr
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np

def mergeResults(fileNameList):
    """
        $BJBNs7W;;$7$?7k2L$r#1$D$N%U%!%$%k$K$^$H$a$k4X?t(B
        fileNameList : $B%U%!%$%kL>$r%j%9%H$K$7$?$b$N!%(B
    """
    df = pd.DataFrame()
    for fileName in fileNameList:
        data = pd.read_csv(fileName)
        df = pd.concat([df, data])
    df = df.reset_index(drop=True)
    return df

def createHeatmap(drugNames, dataLists, subplot, saveName):
    """
        drugNames : csv$B$N%U%!%$%kL>$G;HMQ$7$F$$$kLt:^L>$N%j%9%H(B
        csvdir    : csv$B$,J]B8$5$l$F$$$k%G%#%l%/%H%jL>(B
        subplot   : $B%0%i%U$N?t$r9TNs$G(B[$B2#(B, $B=D(B]
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

        ax.invert_yaxis() # y$B<4$N>e2<$rJQ$($k(B
        setTickLabel(data, ax) # $B<4$N@_Dj(B
        ax.set_ylabel(drugName[0], fontsize = 30) # y$B<4$N%i%Y%k(B
        ax.set_xlabel(drugName[1], fontsize = 30) # x$B<4$N%i%Y%k(B
        ax.tick_params(labelsize=24)

    cbar_ax.tick_params(labelsize = 30)
    fig.tight_layout(rect = [0, 0, .9, 1.])
    plt.savefig(saveName, dpi=300)


def setTickLabel(data, ax):
    a1DoseList = list(set(data["a1"].tolist()))[::20] # y$B<4$KI=<($7$?$$?t$N%j%9%H(B
    a2DoseList = list(set(data["a2"].tolist()))[::20] # x$B<4$KI=<($7$?$$?t$N%j%9%H(B

    # yticks$B$N@_Dj(B
    dataLenY = len(list(set(data["a1"].tolist()))) 
    ax.set_yticks(list(np.linspace(0.5, dataLenY - 0.5, len(a1DoseList))))
    ax.set_yticklabels(list(map(lambda x:"{:.3f}".format(x), a1DoseList)))
    
    # xticks$B$N@_Dj(B
    dataLenX = len(list(set(data["a2"].tolist()))) 
    ax.set_xticks(list(np.linspace(0.5, dataLenX - 0.5, len(a2DoseList))))
    ax.set_xticklabels(list(map(lambda x:"{:.3f}".format(x), a2DoseList)))

if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    drugNameList = itr.combinations_with_replacement(drugNames, 2)
    csvdir = "results/ribo4/csv/sim100"
    
    # merge DataFiles
    for drugName in drugNameList:
        dirName = "{}/{}".format(csvdir, "_".join(drugName))
        fileNameList = ["{}/{}_{}.csv".format(dirName, "_".join(drugName), num) for num in range(101)]
        df = mergeResults(fileNameList)
        df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(drugName)), index=False)

    # SameDrug combination
    drugNameList = [[name, name] for name in drugNames]
    dataLists = ["results/ribo4/csv/sim100/{}_merge.csv".format("_".join(i)) for i in drugNameList]
    doubleSaveName = "results/ribo4/images/sameDrug_sim100.png"
    createHeatmap(drugNameList, dataLists, [2, 2], doubleSaveName)


    # differentDrug combination
    drugNameList = list(itr.combinations(drugNames, 2))
    dataLists = ["results/ribo4/csv/sim100/{}_merge.csv".format("_".join(i)) for i in drugNameList]
    doubleSaveName = "results/ribo4/images/diffDrug_sim100.png"
    createHeatmap(drugNameList, dataLists, [3, 2], doubleSaveName)
