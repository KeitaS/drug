# coding: utf-8
from cynergistic_check import *
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools

dataset = {"dNames": ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"],
           "IC30": {'Streptmycin': 0.9593399047851562, 'Chloramphenicol': 4.3359375, 'Kanamycin': 0.45942234992980957, 'Tetracycline': 0.703125}
           }
a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}

savedir = "images/ribo6"
makedir(savedir)
savedir = "images/ribo6/synergistic"
makedir(savedir)
r_min = 19.3
r_max = 65.8
# IC30を計算
# dataset["IC30"] = calcIC(dataset["dNames"], a_ex, 0.3)


# R30(b,a1,a2) + R50(b) == R30(b^1,a1,a2).R50(b^1)
# R30(a1,b^_free) + a1(b) == R30(a1^1,b).a1(b^1)
# a1 > ~a1

# 単剤
drug_single = [makeDrugDatas("Tetracycline")]
drug_double = [makeDrugDatas("Tetracycline"), makeDrugDatas("Tetracycline")]
drug_double[1]["type"] = "50s"

doses = np.linspace(0, dataset["IC30"]["Tetracycline"] * 2, 100)
result_list = list()
inpData = {"modif": 0}

def single(drug, dose):
    drug[0]["dose"] = dose
    return(run(drug, step=100, inpData={"modif": 0}, legend=["r_u"])[0][-1][1])

def double(drug, dose):
    drug[0]["dose"] = dose
    drug[1]["dose"] = dose
    return(run(drug, step=100, inpData={"modif": 0}, legend=["r_u"])[0][-1][1])

data = pd.DataFrame([[dose, (single(drug_single, dose) - r_min) / r_max, (double(drug_double, dose / 2) - r_min) / r_max] for dose in doses],
                    columns=["dose", "single", "double"])


plt.plot(data["dose"], data["single"], "r-", label="single")
plt.plot(data["dose"], data["double"], "b-", label="double")
plt.legend(loc="upper right")
plt.xlabel("dose")
plt.ylabel("$(r_{u} - r_{min}) / r_{max}$")
plt.savefig("images/ribo6/check_growth/no_lambda_sup.png", dpi=200)
