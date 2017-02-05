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
IC30 = {'Kanamycin': 0.4598379135131836, 'Tetracycline': 0.732421875, 'Chloramphenicol': 4.3212890625, 'Streptmycin': 0.9600967168807983}

# 解離するリアクションがない場合のIC30
# IC30 = {'Tetracycline': 5.2734375, 'Chloramphenicol': 21.09375, 'Streptmycin': 1.46942138671875, 'Kanamycin': 0.6775975227355957}


# 保存するディレクトリの設定
savedir = "results/ribo7"
# makedir(savedir + "/images/test")

# main
# growthPlot(dNames, IC30, 1, savedir+"/images/single_all_growthrate.png")
# heatmap(list(itertools.combinations(dNames, 2)), IC30, 11, savedir + "/images/test/double_heatmap_comb.png", (18, 9))
# heatmap([[dName, dName] for dName in dNames], IC30, 11, savedir + "/images/test/double_heatmap.png", (12, 9))
# neweval(list(itertools.combinations(dNames, 2)), IC30, 11, savedir + "/images/test/double_neweval_comb.png", (18, 9))
# neweval([[dName, dName] for dName in dNames], IC30, 11, savedir + "/images/test/double_neweval.png", (12, 9))
heatmap([[dNames[0], dNames[1]]], IC30, 11, figsize=(12, 9))
