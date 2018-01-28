# coding: utf-8
import matplotlib
matplotlib.use("Agg")
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True
import numpy as np
import pandas as pd
import itertools as itr
import sys
import seaborn as sns
import matplotlib.pylab as plt
from ribo6 import makedir

## liner chart
csvdir = "results/ribo6/csv/liner"
drugNames = ["Streptmycin", "Chloramphenicol"]
imgdir = "results/ribo6/images"
makedir(imgdir)

## normal
for drugName in drugNames:
    normal = pd.read_csv("{}/{}_0.csv".format(csvdir, drugName))
    Lambda_0 = normal["single"][0]
    plt.plot(np.linspace(0, 100, 101), list(map(lambda x: x / Lambda_0, list(normal["single"]))), label="single")
    plt.plot(np.linspace(0, 100, 101), list(map(lambda x: x / Lambda_0, list(normal["double"]))), label="double")
    plt.legend()
    plt.xlabel("Dose", fontsize=20)
    plt.ylabel("Growth Rate", fontsize=20)
    plt.tight_layout()
    plt.savefig("{}/{}_normal.png".format(imgdir, drugName), dpi=300)
    plt.close()

## change
for drugName in drugNames:
    change = pd.read_csv("{}/{}_1.csv".format(csvdir, drugName))
    Lambda_0 = change["single"][0]
    plt.plot(np.linspace(0, 100, 101), list(map(lambda x: x / Lambda_0, list(change["single"]))), label="single")
    plt.plot(np.linspace(0, 100, 101), list(map(lambda x: x / Lambda_0, list(change["double"]))), label="double")
    plt.legend()
    plt.xlabel("Dose", fontsize=20)
    plt.ylabel("Growth Rate", fontsize=20)
    plt.tight_layout()
    plt.savefig("{}/{}_change.png".format(imgdir, drugName), dpi=300)



