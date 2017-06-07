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
# IC30 = calcIC(dNames, {x: 50 for x in dNames}, 0.3)
# print(IC30)

# 薬剤が結合すると解離するリアクションがある場合のIC30
# IC30 = {'Kanamycin': 0.4598379135131836, 'Tetracycline': 0.732421875, 'Chloramphenicol': 4.3212890625, 'Streptmycin': 0.9600967168807983}

# 解離するリアクションがない場合のIC30
IC30 = {'Tetracycline': 5.2734375, 'Chloramphenicol': 21.09375, 'Streptmycin': 1.46942138671875, 'Kanamycin': 0.6775975227355957}


# 保存するディレクトリの設定
csvdir = "results/ribo7/datas/no_reaction/heatmap"
savedir = "results/ribo7/images/no_reaction"
# makedir(csvdir)

# main
# growthPlot(dNames, IC30, 1, savedir+"/images/single_all_growthrate.png")
heatmap(list(itertools.combinations(dNames, 2)), IC30, 11, savedir + "/double_heatmap_comb.png", csvdir) # combination
heatmap([[dName, dName] for dName in dNames], IC30, 11, savedir + "/double_heatmap.png", csvdir) # samedrug
# neweval(list(itertools.combinations(dNames, 2)), IC30, 11, savedir + "/double_neweval_comb.png", (18, 9), csvdir)
# neweval([[dName, dName] for dName in dNames], IC30, 11, savedir + "/double_neweval.png", (12, 9), csvdir)



"""
    ribo5を使って，同じ薬剤を別のターゲットで投与したシミュレーション
"""
# # IC30を計算
# import ribo5
# # IC30 = ribo5.calcIC(dNames, {x: 50 for x in dNames}, 0.3)
# IC30 = {'Tetracycline': 5.2734375, 'Kanamycin': 0.6761401891708374, 'Streptmycin': 1.4652013778686523, 'Chloramphenicol': 21.09375}
# ribo5_diffTargetCheck(dNames, IC30, (12, 9), savedir + "/ribo5_sameDrug_diffTarget.png")
