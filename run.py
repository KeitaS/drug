# coding: utf-8
from ribo7 import *
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools
import math
from modules import *



dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
IC30 = calcIC(dNames, {x: 50 for x in dNames}, 0.3)
print(IC30)

# 薬剤が結合すると解離するリアクションがある場合のIC30
# IC30 = {'Kanamycin': 0.4598379135131836, 'Tetracycline': 0.732421875, 'Chloramphenicol': 4.3212890625, 'Streptmycin': 0.9600967168807983}

# 解離するリアクションがない場合のIC30
IC30 = {'Tetracycline': 5.2734375, 'Chloramphenicol': 21.09375, 'Streptmycin': 1.46942138671875, 'Kanamycin': 0.6775975227355957}


# 保存するディレクトリの設定
# csvdir = "results/ribo5/datas/sharp_flat/original"
# savedir = "results/ribo5/original"
# makedir(csvdir)

# main
# growthPlot(dNames, IC30, 1, savedir+"/images/single_all_growthrate.png")
# heatmap(list(itertools.combinations(dNames, 2)), IC30, 11, savedir + "/double_heatmap_comb.png", csvdir) # combination
# heatmap([[dName, dName] for dName in dNames], IC30, 11, savedir + "/double_heatmap.png", csvdir) # samedrug
# neweval(list(itertools.combinations(dNames, 2)), IC30, 11, savedir + "/double_neweval_comb.png", (18, 9), csvdir)
# neweval([[dName, dName] for dName in dNames], IC30, 11, savedir + "/double_neweval.png", (12, 9), csvdir)

"""
    ribo5を使って，同じ薬剤を別のターゲットで投与したシミュレーション
"""
# IC30を計算
import ribo5

## IC30を計算
# IC30_ribo5 = ribo5.calcIC(dNames, {x: 50 for x in dNames}, 0.3)
IC30_ribo5 = {'Streptmycin': 1.4652013778686523, 'Kanamycin': 0.6761401891708374, 'Chloramphenicol': 21.09375, 'Tetracycline': 5.2734375}

## SterptmycinとChloramphenicolで確認する
# drugNames = list(itertools.combinations_with_replacement(["Streptmycin", "Chloramphenicol"], 2))
#
# plt.figure(figsize=(15, 5))
# num = 1
# target = ["30s", "50s"]
# for drugNameList in drugNames:
#     data = ribo5_targetChange(drugNameList, IC30_ribo5, target=target, modif=0) # ribo5を使って，データ取得
#     data.to_csv("{}/{}_{}_{}_{}.csv".format(csvdir, drugNameList[0], target[0], drugNameList[1], target[1]), index=False) # csvデータ作成
#     plt.subplot(1, 3, num)
#     createHeatmap(data, target)
#     num += 1
# plt.tight_layout()
# plt.savefig("{}/sharp_flat.png".format(savedir), dpi=300)


"""
    カーブが急なものとカーブがなだらかなもので比較．
"""
# csvdir = "results/ribo7/datas/default/sharp_flat"
# savedir = "results/ribo7/images/default"
#
#
# drugNames = list(itertools.combinations_with_replacement(["Streptmycin", "Chloramphenicol"], 2))
# plt.figure(figsize=(15, 10))
# num = 1
#
# for drugNameList in drugNames:
#     drugs = [makeDrugDatas(drugNameList[0]), makeDrugDatas(drugNameList[1])]
#     drugs[0]["target"] = "r30"
#     drugs[1]["target"] = "r30"
#     doses = [[x, y] for x in np.linspace(0, IC30[drugNameList[0]] * 2, 11) for y in np.linspace(0, IC30[drugNameList[1]] * 2, 11)]
#     data = pd.DataFrame([[round(dose[0], 2), round(dose[1], 2), sim(drugs, dose)] for dose in doses],
#                         columns=["a1", "a2", "growth"]) # simulation & create data
#     data.to_csv("{}/{}_{}_{}_{}.csv".format(csvdir, drugs[0]["name"], drugs[0]["target"], drugs[1]["name"], drugs[1]["target"]), index=False)
#     plt.subplot(2, 3, num)
#     createHeatmap(data, ["30s", "30s"])
#     num += 1
#
# for drugNameList in drugNames:
#     drugs = [makeDrugDatas(drugNameList[0]), makeDrugDatas(drugNameList[1])]
#     drugs[0]["target"] = "r30"
#     drugs[1]["target"] = "r50"
#     doses = [[x, y] for x in np.linspace(0, IC30[drugNameList[0]] * 2, 11) for y in np.linspace(0, IC30[drugNameList[1]] * 2, 11)]
#     data = pd.DataFrame([[round(dose[0], 2), round(dose[1], 2), sim(drugs, dose)] for dose in doses],
#                         columns=["a1", "a2", "growth"]) # simulation & create data
#     data.to_csv("{}/{}_{}_{}_{}.csv".format(csvdir, drugs[0]["name"], drugs[0]["target"], drugs[1]["name"], drugs[1]["target"]), index=False)
#     plt.subplot(2, 3, num)
#     createHeatmap(data, ["30s", "50s"])
#     num += 1
#
# plt.tight_layout()
# plt.savefig("{}/sharp_flat.png".format(savedir), dpi=300)

"""
ribo5を用いて，モデル変更前のデータを出力
"""
csvdir = "results/ribo5/datas/original"
savedir = "results/ribo5/images/original"
makedir(csvdir)

# ribo5_heatmap(list(itertools.combinations(dNames, 2)), IC30_ribo5, savename="{}/double_heatmap_comb.png".format(savedir), csvdir=csvdir)
# ribo5_heatmap([[dName, dName] for dName in dNames], IC30_ribo5, savename="{}/double_heatmap.png".format(savedir), csvdir=csvdir)
