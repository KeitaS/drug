# coding: utf-8
from ribo5 import *
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools

savedir = "images/ribo5"
makedir(savedir)

dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}

# IC30 = calcIC(dNames, a_ex, .3) # IC30を計算
IC30 = {'Kanamycin': 0.6761398315429688, 'Streptmycin': 1.4652465820312497, 'Chloramphenicol': 22.5, 'Tetracycline': 5.25}
modif2_K_ma = {'Kanamycin': 8.9, 'Streptmycin': 9.2, 'Chloramphenicol': 23.1, 'Tetracycline': 24.7}

drug_comb = list(itertools.combinations(dNames, 2))

# create Growth heatmap
def createGrowthHeatmap(modif):
    plt.figure(figsize=(12, 6))
    drug_comb = [[i, i] for i in dNames]
    for index, dName in enumerate(drug_comb):
        drugs = [makeDrugDatas(dName[0]), makeDrugDatas(dName[1])]
        doses = [np.linspace(0, IC30[dName[0]] * 2, 11), np.linspace(0, IC30[dName[1]] * 2, 11)]
        doses = list(itertools.product(doses[0], doses[1]))
        result_list = []
        inpData = {"modif": modif}

        for dose in doses:
            result = doseResponse(drugs, dose, inpData=inpData)
            result_list.append([round(dose[0], 2), round(dose[1], 2), result])

        data = pd.DataFrame(result_list, columns=["a1", "a2", "growth"])
        # data = pd.DataFrame(result_list, columns=[dName[0], dName[1], "growth"])
        plt.subplot(2, 3, index + 1)
        growthHeatmap(data=data, values="growth", index="a1", columns="a2", title=dName[0])

    plt.tight_layout()
    plt.savefig("{}/modif2/modification_2_double.png".format(savedir), dpi=200)
    # plt.show()
    plt.close()


# create Epsilon heatmap
def createEpsilonHeatmap(modif):
    plt.figure(figsize=(12, 6))
    cmap = makeCmap()
    slope_list = [1./4, 1./2, 1., 2., 4.] # 傾き

    for index, dName in enumerate(drug_comb):
        midPointList = midPointCombination([IC30[dName[0]] * 2, IC30[dName[1]] * 2])
        data = []
        linertype = 0
        inpData = {"modif": modif}
        for slope in slope_list:
            for pCount, midPoint in enumerate(midPointList):
                doses = [midPoint[0] * (1 + slope), midPoint[1] * (1 + (1 / slope))]
                linertype = calcEpsilon(dName, doses, inpData=inpData)
                data.append([slope, (pCount + 1) * 10, linertype])

        plt.subplot(2, 3, index + 1)
        data = pd.DataFrame(data, columns=["S", "I", "growth_type"])
        evalHeatmap(data, cmap, "growth_type", "I", "S", ylabel="MidPoint", xlabel="Slope", title="{} vs {}".format(dName[0][:3], dName[1][:3]))

    plt.tight_layout()
    plt.savefig("{}/modification_comb_{}_oldeval.png".format(savedir, modif), dpi=200)
    # plt.show()
    plt.close()

# create neweval heatmap
def createNewevalHeatmap(modif, norm):
    plt.figure(figsize=(12, 6))
    cmap = generate_cmap(["mediumblue", "white", "orangered"])
    slope_list = [1./4, 1./2, 1., 2., 4.] # 傾き\
    drug_comb = [[i, i] for i in dNames]

    for index, dName in enumerate(drug_comb):
        drugs = [makeDrugDatas(dName[0]), makeDrugDatas(dName[1])]
        midPointList = midPointCombination([IC30[dName[0]] * 2, IC30[dName[1]] * 2])
        data = []
        linertype = 0
        inpData = {"modif": modif}
        for slope in slope_list:
            for pCount, midPoint in enumerate(midPointList):
                result_list = []
                doses = createSlopedose(slope, midPoint)
                for dose in doses:
                    result_list.append(doseResponse(drugs, dose, inpData=inpData))

                if norm == True:
                    buffpoint = calcBufferingPoint([dName[0], dName[1]], [doses[-1][0], doses[0][1]]) # buffering point を計算
                    linertype = checkLinerType(result_list, 1.0e-6, 2, buffpoint)
                else:
                    linertype = checkLinerType(result_list, 1.0e-6, 1)
                data.append([slope, (pCount + 1) * 10, linertype])

        plt.subplot(2, 3, index + 1)
        data = pd.DataFrame(data, columns=["S", "I", "growth_type"])
        evalHeatmap(data, cmap, "growth_type", "I", "S", ylabel="MidPoint", xlabel="Slope", title=dName)
        # evalHeatmap(data, cmap, "growth_type", "I", "S", ylabel="MidPoint", xlabel="Slope", title="{} vs {}".format(dName[0][:3], dName[1][:3]))

    plt.tight_layout()
    if norm == True: plt.savefig("{}/modif2/modification_{}_neweval_norm.png".format(savedir, modif), dpi=200)
    else: plt.savefig("{}/modif2/modification_{}_neweval.png".format(savedir, modif), dpi=200)
    # if norm == True: plt.savefig("{}/modification_comb_{}_neweval_norm.png".format(savedir, modif), dpi=200)
    # else: plt.savefig("{}/modification_comb_{}_neweval.png".format(savedir, modif), dpi=200)
    # plt.show()
    plt.close()

# 実行
createNewevalHeatmap(2, True)
createNewevalHeatmap(2, False)
