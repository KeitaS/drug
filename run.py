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
modif2_K_ma = {'Kanamycin': 8.9, 'Streptmycin': 9.2, 'Chloramphenicol': 23.1, 'Tetracycline': 24.7}


# create Growth heatmap
plt.figure(figsize=(12, 9))

for index, dName in enumerate(dNames):
    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
    doses = np.linspace(0, IC30[dName] * 2, 11)
    doses = list(itertools.product(doses, doses))
    result_list = []
    inpData = {"K_ma": modif2_K_ma[dName], "modif": 2}

    for dose in doses:
        result = doseResponse(drugs, dose, inpData=inpData)
        result_list.append([round(dose[0], 2), round(dose[1], 2), result])

    data = pd.DataFrame(result_list, columns=["a1", "a2", "growth"])
    plt.subplot(2, 2, index + 1)
    growthHeatmap(data=data, values="growth", index="a1", columns="a2", title=dName)

plt.tight_layout()
plt.savefig("{}/modification_2_double.png".format(savedir), dpi=200)
# plt.show()
plt.close()


# create Epsilon heatmap
"""
plt.figure(figsize=(12, 9))
cmap = makeCmap()
slope_list = [1./4, 1./2, 1., 2., 4.] # 傾き

for index, dName in enumerate(dNames):
    doses = np.linspace(0, IC30[dName] * 2, 11)
    midPointList = [[doses[i] / 2, doses[i] / 2] for i in range(len(doses)) if i > 0] # 中点のリスト
    data = []
    linertype = 0
    for slope in slope_list:
        for pCount, midPoint in enumerate(midPointList):
            doses = [midPoint[0] * (1 + slope), midPoint[1] * (1 + (1 / slope))]
            linertype = calcEpsilon([dName, dName], doses, inpData={"K_ma": modif2_K_ma[dName], "modif": 2})
            data.append([slope, (pCount + 1) * 10, linertype])

    plt.subplot(2, 2, index + 1)
    data = pd.DataFrame(data, columns=["S", "I", "growth_type"])
    evalHeatmap(data, cmap, "growth_type", "I", "S", ylabel="MidPoint", xlabel="Slope", title=dName)

plt.tight_layout()
plt.savefig("{}/modification_2_oldeval.png".format(savedir), dpi=200)
# plt.show()
plt.close()
"""

# create neweval heatmap
"""
plt.figure(figsize=(12, 9))
cmap = generate_cmap(["mediumblue", "white", "orangered"])
slope_list = [1./4, 1./2, 1., 2., 4.] # 傾き


for index, dName in enumerate(dNames):
    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
    doses = np.linspace(0, IC30[dName] * 2., 4)
    midPointList = [[doses[i] / 2., doses[i] / 2.] for i in range(len(doses)) if i > 0] # 中点のリスト
    data = []
    linertype = 0
    for slope in slope_list:
        for pCount, midPoint in enumerate(midPointList):
            result_list = []
            doses = createSlopedose(slope, midPoint)
            for dose in doses:
                result_list.append(doseResponse(drugs, dose, inpData={"K_ma": modif2_K_ma[dName], "modif": 2}))

            # buffpoint = calcBufferingPoint([dName, dName], [doseX[-1], doseY[0]]) # buffering point を計算
            # linertype = checkLinerType(result_list, 1.0e-6, 2, buffpoint)
            linertype = checkLinerType(result_list, 1.0e-6, 1)
            data.append([slope, (pCount + 1) * 10, linertype])

    plt.subplot(2, 2, index + 1)
    data = pd.DataFrame(data, columns=["S", "I", "growth_type"])
    print data
    evalHeatmap(data, cmap, "growth_type", "I", "S", ylabel="MidPoint", xlabel="Slope", title=dName)

plt.tight_layout()
# plt.savefig("{}/modification_2_neweval.png".format(savedir), dpi=200)
plt.show()
plt.close()
"""

# ①に近くなるK_maを特定する
"""
## ①のときのgrowthを計算
for index, dName in enumerate(dNames):
    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
    doses = [np.linspace(0, IC30[dName] * 2, 11), np.linspace(IC30[dName] * 2, 0, 11)]

modif1 = {}
for index, dName in enumerate(dNames):
    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
    K_ma = 15.
    point = np.linspace(0, IC30[dName] * 2, 11)[7]
    doseX = np.linspace(0, point, 11)
    doseY = np.linspace(point, 0, 11)
    doses = [[doseX[i], doseY[i]] for i in range(len(doseX))]
    result_list = []
    for dose in doses:
        drugs[0]["dose"] = dose[0]
        drugs[1]["dose"] = dose[1]
        result, legend = run(drugs, step=100, inpData={"K_ma": K_ma, "modif": 1}, legend=["r_u"])
        result = calcGrowthrate(result[-1][1])
        result_list.append(result)
    modif1[dName] = max(result_list)

print "modif1 Complete"

## ②を実行して，①と比較
K_ma_list = {}
for index, dName in enumerate(dNames):
    print dName
    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
    K_ma = 15.
    point = np.linspace(0, IC30[dName] * 2, 11)[7] # 一番顕著に出てくる部分で比較
    count = 0
    check_flag = 0
    while K_ma < 100 and K_ma > 0:
        print K_ma

        doseX = np.linspace(0, point, 11)
        doseY = np.linspace(point, 0, 11)
        doses = [[doseX[i], doseY[i]] for i in range(len(doseX))]
        result_list = []
        for dose in doses:
            drugs[0]["dose"] = dose[0]
            drugs[1]["dose"] = dose[1]
            result, legend = run(drugs, step=100, inpData={"K_ma": K_ma, "modif": 2}, legend=["r_u"])
            result = calcGrowthrate(result[-1][1])
            result_list.append(result)
            if result - modif1[dName] > 0.001: # modif1より0.01より大きくなったらbreak
                break
        print modif1[dName]
        print max(result_list)

        if abs(modif1[dName] - max(result_list)) < 0.001: # modif1に近かったら
            break
        elif modif1[dName] > max(result_list): # modif1より小さかったら
            K_ma -= 0.1
            if check_flag == -1: count = 0
            else:
                check_flag = -1
                count += 1

        else:
            K_ma += 0.1
            if check_flag == 1: count = 0
            else:
                check_flag = 1
                count += 1

        if count > 5:
            break

    modif2_K_ma[dName] = round(K_ma, 1)
print modif2_K_ma

with open("output.txt", "w") as wf:
    for i, j in K_ma_list.items():
        wf.write("{}, {}\n".format(i, j))
## LinerChartを作成

for index, dName in enumerate(dNames):
    drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
    point = np.linspace(0, IC30[dName] * 2, 11)[7] / 2
    doses = createSlopedose(1, [point, point])
    data = []
    for count, dose in enumerate(doses):
        data.append([count, doseResponse(drugs, dose, inpData={"K_ma": 15., "modif": 1})])
        data[count].append(doseResponse(drugs, dose, inpData={"K_ma": modif2_K_ma[dName], "modif": 2}))
    data = np.array(data)
    plt.subplot(2, 2, index + 1)
    plt.plot(data.T[0], data.T[1], label="modif1(K_ma=15)")
    plt.plot(data.T[0], data.T[2], label="modif2(K_ma={})".format(modif2_K_ma[dName]))
    plt.ylabel("growth")
    plt.xlabel("count")
    plt.title(dName)
    plt.legend(loc="upper right")

plt.tight_layout()
plt.savefig("{}/check_modif1_and_modif2_growth.png".format(savedir), dpi=200)
# plt.show()
plt.close()
"""

# nonlinear vs nonlinear combination
"""
dName = "Streptmycin"
drugs = [makeDrugDatas(dName), makeDrugDatas(dName)]
drugs[1]["type"] = "50s"
doses = np.linspace(0, IC30[dName] * 2, 11)
doses = list(itertools.product(doses, doses))
result_list = []
for dose in doses:
    result_list.append([round(dose[0], 2), round(dose[1], 2), doseResponse(drugs, dose, inpData={"modif": 0})])
data = pd.DataFrame(result_list, columns=["30s", "50s", "growth"])
heatmap = pd.pivot_table(data=data, values="growth", index="30s", columns="50s")
sns.heatmap(heatmap)
plt.title("30s vs 50s(non linear)")
plt.ylabel("30s")
plt.xlabel("50s")
plt.tick_params(labelsize=7)
plt.savefig("{}/heatmap_nonlinear_difftype.png".format(savedir), dpi=200)
# plt.show()
plt.close()
"""
