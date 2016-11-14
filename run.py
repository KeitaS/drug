# coding: utf-8
from ribo5 import *
from modules import *
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt

savedir = "images/ribo5"
makedir(savedir)

a_ex = 0.7 * 2

# create line chart
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
