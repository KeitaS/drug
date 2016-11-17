# coding: utf-8
from ribo5 import *
from modules import *
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

# check in linechart
"""
## single
drugs = [makeDrugDatas("Kanamycin")]
doses = np.linspace(0, a_ex, 31)
result_list = []
for dose in doses:
    drugs[0]["dose"] = dose
    result, legend = run(drugs, step=100, legend=["r_u"])
    result = calcGrowthrate(result[-1][1])
    result_list.append([dose, result])

## double
drugs = [makeDrugDatas("Kanamycin"), makeDrugDatas("Kanamycin")]
doses = np.linspace(0, a_ex / 2, 31)
for index, dose in enumerate(doses):
    drugs[0]["dose"] = dose
    drugs[1]["dose"] = dose
    result, legend = run(drugs, step=100, legend=["r_u"])
    result = calcGrowthrate(result[-1][1])
    result_list[index].append(result)
result_list = np.array(result_list)
plt.plot(result_list.T[0], result_list.T[1], label="single")
plt.plot(result_list.T[0], result_list.T[2], label="double")
plt.legend()
plt.show()
"""

# create heatmap
"""
plt.figure(figsize=(6, 9))
for index, dName in enumerate(dNames):
    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
    doses = np.linspace(0, IC30[dName] * 2, 11)
    doses = list(itertools.product(doses, doses))
    result_list = []

    for dose in doses:
        drugs[0]["dose"] = dose[0]
        drugs[1]["dose"] = dose[1]
        result, legend = run(drugs, step=100, inpData={"K_ma": 30.}, legend=["r_u"])
        result = calcGrowthrate(result[-1][1])
        result_list.append([round(dose[0], 2), round(dose[1], 2), result])
    data = pd.DataFrame(result_list, columns=["a1", "a2", "growth"])
    heatmap = pd.pivot_table(data=data, values="growth", index="a1", columns="a2")
    plt.subplot(2,2,index + 1)
    sns.heatmap(heatmap)
    plt.title(dName)
    plt.ylabel("a1")
    plt.xlabel("a2")
    plt.tick_params(labelsize=7)
plt.tight_layout()
# plt.savefig("{}/modification_1_double_10.png".format(savedir), dpi=200)
plt.show()
plt.close()
"""


# 指標のヒートマップ
plt.figure(figsize=(12, 9))
cmap = generate_cmap(["mediumblue", "white", "orangered"])
slope_list = [1./4, 1./2, 1, 2, 4] # 傾き

for index, dName in enumerate(dNames):
    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
    doses = np.linspace(0, IC30[dName] * 2, 11)
    midPointList = [[doses[i] / 2, doses[i] / 2] for i in range(len(doses)) if i > 0] # 中点のリスト
    data = []
    K_ma = 15.
    linertype = 0
    for slope in slope_list:
        for pCount, midPoint in enumerate(midPointList):
            result_list = []
            doseX = np.linspace(0, midPoint[0] * (1 + slope), 11)
            doseY = np.linspace(0, midPoint[1] * (1 + (1 / slope)), 11)[::-1]
            doses = [[doseX[i], doseY[i]] for i in range(len(doseX))]
            for dose in doses:
                drugs[0]["dose"] = dose[0]
                drugs[1]["dose"] = dose[1]
                result, legend = run(drugs, step=100, inpData={"K_ma": K_ma}, legend=["r_u"])
                result_list.append(calcGrowthrate(result[-1][1]))

            buffpoint = calcBufferingPoint([dName, dName], [doseX[-1], doseY[0]]) # buffering point を計算
            linertype = checkLinerType(result_list, 1.0e-6, 2, buffpoint)
            data.append([slope, (pCount + 1) * 10, linertype])

    data = pd.DataFrame(data, columns=["S", "I", "growth_type"])

    heatmap = pd.pivot_table(data=data, values="growth_type", index="I", columns="S") # x軸を0, y軸を1番目の薬剤にしてグラフデータ化
    plt.subplot(2, 2, index + 1) # 1つの画像データに集約
    sns.heatmap(heatmap, vmin=-1, vmax=1, cmap=cmap, linewidths=.3, annot=True, annot_kws={"size": 7})
    ax = plt.gca()
    ax.invert_yaxis() # ヒートマップのy軸の逆転
    plt.tick_params(labelsize=7)
    plt.ylabel("MidPoint")
    plt.xlabel("Slope")
    plt.title(dName)

plt.tight_layout()
plt.savefig("{}/modification_1_antagocheck_15.png".format(savedir), dpi=200)
# plt.show()
plt.close()


# oldeval heatmap
"""
plt.figure(figsize=(12, 9))
cmap = makeCmap()
slope_list = [1./4, 1./2, 1, 2, 4] # 傾き

for index, dName in enumerate(dNames):
    doses = np.linspace(0, IC30[dName] * 2, 11)
    midPointList = [[doses[i] / 2, doses[i] / 2] for i in range(len(doses)) if i > 0] # 中点のリスト
    data = []
    K_ma = 15.
    linertype = 0
    for slope in slope_list:
        for pCount, midPoint in enumerate(midPointList):
            ## 単剤
            doses = [midPoint[0] * (1 + slope), midPoint[1] * (1 + (1 / slope))]
            result_list = []
            for i in range(3):
                if i < 2:
                    drugs = [makeDrugDatas(dName)]
                    drugs[0]["dose"] = doses[i]
                    result, legend = run(drugs, step=100, inpData={"K_ma": K_ma}, legend=["r_u"])
                    result = calcGrowthrate(result[-1][1])
                    result_list.append(result)
                else:
                    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
                    drugs[0]["dose"] = doses[0]
                    drugs[1]["dose"] = doses[1]
                    result, legend = run(drugs, step=100, inpData={"K_ma": K_ma}, legend=["r_u"])
                    result = calcGrowthrate(result[-1][1])
                    result_list.append(result)
            linertype = epsilon(result_list[0], result_list[1], result_list[2])
            data.append([slope, (pCount + 1) * 10, linertype])

    data = pd.DataFrame(data, columns=["S", "I", "growth_type"])

    heatmap = pd.pivot_table(data=data, values="growth_type", index="I", columns="S") # x軸を0, y軸を1番目の薬剤にしてグラフデータ化
    plt.subplot(2, 2, index + 1) # 1つの画像データに集約
    sns.heatmap(heatmap, vmin=-1, vmax=1, cmap=cmap, linewidths=.3, annot=True, annot_kws={"size": 7})
    ax = plt.gca()
    ax.invert_yaxis() # ヒートマップのy軸の逆転
    plt.tick_params(labelsize=7)
    plt.ylabel("MidPoint")
    plt.xlabel("Slope")
    plt.title(dName)

plt.tight_layout()
plt.savefig("{}/modification_1_oldeval_15.png".format(savedir), dpi=200)
plt.close()
"""
