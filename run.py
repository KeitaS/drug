# coding: utf-8
from ribo6 import *
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools

dataset = {"dNames": ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"],
           "IC30": {'Kanamycin': 0.6761398315429688, 'Streptmycin': 1.4652465820312497, 'Chloramphenicol': 22.5, 'Tetracycline': 5.25}
           }
a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}

savedir = "images/ribo6"
makedir(savedir)
savedir = "images/ribo6/synergistic"
makedir(savedir)

# R30(b,a1,a2) + R50(b) == R30(b^1,a1,a2).R50(b^1)
# R30(a1,b^_free) + a1(b) == R30(a1^1,b).a1(b^1)
# a1 > ~a1

createGrowthHeatmap(dataset=dataset, modif=0, savename=savedir+"/heatmap.png", comb=False)
print("growth no comb.")
createEpsilonHeatmap(dataset=dataset, modif=0, savename=savedir+"/oldeval.png", comb=False)
print("epsilon no comb.")
createNewevalHeatmap(dataset=dataset, modif=0, norm=False, savename=savedir+"/neweval.png", comb=False)
print("neweval no comb.")
createNewevalHeatmap(dataset=dataset, modif=0, norm=True, savename=savedir+"/neweval_norm.png", comb=False)
print("neweval no comb normalize.")
createGrowthHeatmap(dataset=dataset, modif=0, savename=savedir+"/heatmap_comb.png", comb=True)
print("growth comb.")
createEpsilonHeatmap(dataset=dataset, modif=0, savename=savedir+"/oldeval_comb.png", comb=True)
print("epsilon comb.")
createNewevalHeatmap(dataset=dataset, modif=0, norm=False, savename=savedir+"/neweval_comb.png", comb=True)
print("neweval comb.")
createNewevalHeatmap(dataset=dataset, modif=0, norm=True, savename=savedir+"/neweval_norm_comb.png", comb=True)
print("neweval comb normalize.")
