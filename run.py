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
plt.figure(figsize=(12, 9))
for index, dName in enumerate(dNames):
    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
    doses = np.linspace(0, IC30[dName] * 2, 11)
    doses = list(itertools.product(doses, doses))
    result_list = []

    for dose in doses:
        drugs[0]["dose"] = dose[0]
        drugs[1]["dose"] = dose[1]
        result, legend = run(drugs, step=100, inpData={"K_ma": 1.}, legend=["r_u"])
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
plt.savefig("{}/modification_1_double_1.png".format(savedir), dpi=200)
plt.close()
