#coding: utf-8

import seaborn as sns
import pandas as pd
import matplotlib.pylab as plt
import itertools as itr

if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    dataPath = "results/ribo8/single/csv"

    plt.figure(figsize=(16, 12))
    for index, drugName in enumerate(drugNames):
        data = pd.read_csv("{}/{}.csv".format(dataPath, drugName))
        plt.subplot(2, 2, index + 1)
        plt.plot(data["dose"], data["growth"])
        plt.xlabel("Dose", fontsize=20)
        plt.ylabel("Growth Rate", fontsize=20)
        plt.title(drugName, fontsize=30)
        plt.tick_params(labelsize=14)
    plt.tight_layout()
    plt.savefig("results/ribo8/single/image/result.png", dpi=300)
