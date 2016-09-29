#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True
from ribo4 import *

# 保存用ディレクトリの作成
savedir = "./images/result4"
makedir(savedir)

# 初期値
r_min = 19.3 # Lambdaの計算で使用
K_t = 6.1 * 10 ** -2 # Lambdaの計算で使用
Kd = 6. # r30とr50の解離率
p = .1 # r30とr50の比率

## drug data
dNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
Lambda_0 =  [1.35, 0.85, 0.40] # 1/h (1.35, 0.85, 0.40)
Lambda_0_a = {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83} # 1/h
IC50 = {"Streptmycin": [0.41, 0.28, 0.196], "Kanamycin": [0.246, 0.096, 0.065], "Tetracycline": [0.5, 0.6, 1.45], "Chloramphenicol": [2.85, 2.65, 5.7]} # µg/ml
IC50_a = {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49} # µg/ml
a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}
dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}
