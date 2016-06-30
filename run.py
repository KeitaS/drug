# coding: utf-8

"""
"""

from ribo import run as run1
from ribo3 import run as run3
from ribo3 import makeGraph as mg
import numpy as np

# 保存用ディレクトリの作成
import os
if not os.path.exists("images/result"):
    os.makedirs("images/result")
del(os)

r_min = 19.3
K_t = 6.1 * 10 ** -2
Kd = 5.
p = 0.1

## drug data
drug = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
Lambda_0 =  [1.35, 0.85, 0.40] # 1/h (1.35, 0.85, 0.40)
Lambda_0_a = {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83} # 1/h
IC50 = {"Streptmycin": [0.41, 0.28, 0.196], "Kanamycin": [0.246, 0.096, 0.065], "Tetracycline": [0.5, 0.6, 1.45], "Chloramphenicol": [2.85, 2.65, 5.7]} # µg/ml
IC50_a = {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49} # µg/ml
a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}

## Kanamycinで比較
name = "Kanamycin"
# dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}
legend = ["r_u"]
N_IC50 = {"30s": 0, "50s": 0, "ribo": 0}


# 今までのモデルと、現在のモデルの比較
## ribo.py
dataset = {"Lambda_0": Lambda_0[0], "Lambda_0_a": Lambda_0_a[name], "IC50": IC50[name][0], "IC50_a":IC50_a[name]}

all_result = []
for dose in np.linspace(0, a_ex[name], 51):
    result = run1(dose, dataset=dataset)
    all_result.append([dose, result["growth"] / Lambda_0[0]])

## ribo3.py
dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}

drugs = [{"name": name, "type": "ribo", "dose": .0, "Lambda_0_a": Lambda_0_a[name], "IC50": IC50[name][0], "IC50_a": IC50_a[name]}]
drug_types = ["30s", "50s", "ribo"]

for num in range(len(drug_types)):
    drug_type = drug_types[num]
    drugs[0]["type"] = drug_type
    count = 0

    # doseを振り、モデルを実行
    for dose in np.linspace(0, a_ex[name], 51):
        print dose
        drugs[0]["dose"] = dose
        result, legend = run3(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0] # 結果をgrowthに書き換え
        all_result[count].append(result)
        count += 1

    # グラフ描画
    legend = ["ribo1", "30s", "50s", "ribo"]
    xlabel = "Extracellular Antibiotic Concentration $a_{ex}$($\mu g/ml$)"
    ylabel = "Normalized Growth Rate($\lambda/\lambda_{0}$)"
    title = "Single Dose Response of Kanamycin"

mg(np.array(all_result), savename="images/ribo_ribo3.png", legend=legend, title=title, xlabel=xlabel, ylabel=ylabel)
