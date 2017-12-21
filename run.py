#coding: utf-8

import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools as itr

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
    plt.savefig("results/ribo8/single/image/result.png", dpi=300)

def createHeatmap(data, drugNames, cbar=False, cmap=False):
    if not cmap: cmap = sns.diverging_palette(220, 10, as_cmap=True) # coler map
    heatmapData = pd.pivot_table(data=data, values="growth", index="dose1", columns="dose2") # heatmap data
    ax = sns.heatmap(heatmapData, cbar=cbar, cmap=cmap, square=True) # create Heatmap
    ax.invert_yaxis()

    setTickLabel(data, ax)

    ax.set_ylabel(drugNames[0], fontsize=16) # create ylabel
    ax.set_xlabel(drugNames[1], fontsize=16) # create xlabel

    ## virtual drug
    # if drugNames[0] == "Streptmycin": ax.set_ylabel("Pattern A", fontsize=axizFontSize) # create ylabel
    # else : ax.set_ylabel("Pattern B", fontsize=axizFontSize) # create ylabe
    # if drugNames    [1] == "Streptmycin": ax.set_xlabel("Pattern A", fontsize=axizFontSize) # create xlabel
    # else: ax.set_xlabel("Pattern B", fontsize=axizFontSize) # create xlabel

    ax.tick_params(labelsize=16)


def setTickLabel(data, ax):
    a1DoseList = list(set(data["dose1"].tolist()))[::20] # y$B<4$KI=<($7$?$$?t$N%j%9%H(B
    a2DoseList = list(set(data["dose2"].tolist()))[::20] # x$B<4$KI=<($7$?$$?t$N%j%9%H(B

    # yticks$B$N@_Dj(B
    dataLenY = len(list(set(data["dose1"].tolist()))) 
    ax.set_yticks(list(np.linspace(0.5, dataLenY - 0.5, len(a1DoseList))))
    ax.set_yticklabels(list(map(str, a1DoseList)))
    
    # xticks$B$N@_Dj(B
    dataLenX = len(list(set(data["dose2"].tolist()))) 
    ax.set_xticks(list(np.linspace(0.5, dataLenX - 0.5, len(a2DoseList))))
    ax.set_xticklabels(list(map(str, a2DoseList)))


if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    
    # single simulation 
    singleDataPath = "results/ribo8/single/csv"
    singleSaveName = "results/ribo8/single/image"

    # combinatorial simulation
    ## merge DataFiles
    drugNameList = itr.combinations_with_replacement(drugNames, 2)
    for drugName in drugNameList:
        fileNameList = ["results/ribo8/double/normal/{}_{}.csv".format("_".join(drugName), num) for num in range(101)]
        df = mergeResults(fileNameList)
        df.to_csv("results/ribo8/double/normal/{}_merge.csv".format("_".join(drugName)), index=False)

    ## SameDrug combination
    doubleSaveName = "resutls/ribo8/double/image/sameDrug.png"
    drugNameList = [[name, name] for name in drugNames]

    plt.figure(figsize=(20, 20))
    for index, drugName in enumerate(drugNameList):
        plt.subplot(2, 2, index + 1)
        data = pd.read_csv("results/ribo8/double/normal/{}_merge.csv".format("_".join(drugName)))
        createHeatmap(data, drugNames)
    plt.tight_layout()
    plt.savefig(doubleSaveName, dpi=300)

    ## differentDrug combination
    doubleSaveName = "resutls/ribo8/double/image/diffDrug.png"
    drugNameList = itr.combinations(dNames, 2)
    
    plt.figure(figsize=(20, 30))
    for index, drugName in enumerate(drugNameList):
        plt.subplot(2, 3, index + 1)
        data = pd.read_csv("results/ribo8/double/normal/{}_merge.csv".format("_".join(drugName)))
        createHeatmap(data, drugNames)
    plt.tight_layout()
    plt.savefig(doubleSaveName, dpi=300)
   

