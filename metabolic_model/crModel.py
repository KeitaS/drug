#coding: utf-8

# ribo
from ribo import run

# cobra
from cobra.test import create_test_model
from multi_knockdown import multi_knockdown
from cobra.io import read_sbml_model
import re

# drugbank
from drug2bnum import drug2bnum

# etc
from collections import defaultdict

def checkinCOBRA(genes, model):
    """
    12/21 
    あるものとないものどちらもreturnする
    cobraに入っている遺伝子の検索用モジュール

    genes: 検索したい遺伝子のリスト
    model: COBRAで作成したモデル
    """
    modelgenes = model.genes
    result_true = []
    result_false = []
    for gene in genes:
        if gene in modelgenes:
            result_true.append(gene)
        else:
            result_false.append(gene)
    
    return result_true, result_false


def makeRiboData(model="model/ribo_data.csv"):
    """
    drug idから初期値等のデータを抽出
    return ... {drug_id:{"data":{"Lambda_0_a":, "IC50":, "IC50_a":,}, drug name, b-number}}
    """
    ribodata = {}
    with open(model, "r") as fopen:
        for line in fopen.readlines():
            if re.match("#", line):
                line = line.strip("# ")
                label = re.split(",", line.strip()) # label

            else:
                line = re.split(",", line.strip())
                for index, val in enumerate(line):
                    if index == 0:
                        key = val
                        ribodata[key] = {"data": {}}
                    else:
                        if label[index] == "b-number":
                            val = re.split(";", val)
                            ribodata[key][label[index]] = val
                        
                        elif label[index] != "data": 
                            ribodata[key][label[index]] = val
                        
                        else:
                            ribodata[key]["data"][label[index]] = float(val)
    
    return ribodata



def crModel(drugs, step_int=100, model_c="model/model.xml", model_r="model/ribo_data.csv", drug_data="model/approved_target_ids_all.csv"):
    """
    12/21
    cobraとriboの複合モデル
    COBRAの判別を目指す
    drugsでa_exも渡せるように

    drugs: {DrugBank ID: Dose}
    step_int: riboModel step
    model_c: COBRA model
    model_r: riboModel data
    drug_data: drug data was downloaded from DrugBank. 

    riboのDose調整にはstreptmycin

    ToDo:
        riboのDose調整を全薬剤で適応
        riboのinputを詳細化
        riboからcobraのフィードバック
    """

    # load phase
    model_c = read_sbml_model(model_c)
    model_c.optimize(solver="glpk")
    wt_flux = model_c.solution.x_dict.copy()
    wt_model = model_c.copy()
    drugs = {drug:{"a_ex": drugs[drug]} for drug in drugs.keys()} # drug名をkeyにしたdictになる
    ribo_bnum = ["b3230", "b3319", "b3342", "b3296", "b3298", "P61177", "b3307", "P60725", "b3316", "b3313"] # riboをターゲットとしているbnumber, ここを増やしていく.

    # check phase
    drug_data = drug2bnum(drug_data)
    dataset_r = makeRiboData(model_r) # {drug_id:{"data":{"Lambda_0_a":, "IC50":, "IC50_a":,}, drug name, b-number}}
    
    cobra_target = defaultdict(lambda: 0.0) # cobraをターゲットにする遺伝子とfold_changeのdict
    ribo_target = [] # riboをターゲットにする薬剤のリスト

    # リボソームは単一薬剤のみなので、ターゲットになっていたらfragをTrueにする
    ribo_data = {"flag": False, "a_ex": 0, "dataset": {}}
    
    for drug in drugs.keys():
        if drug_data.get(drug):
            drugs[drug]["all"] = drug_data[drug]
        else:
            drugs[drug]["all"] = None

        drugs[drug]["cobra"] = None
        drugs[drug]["ribo"] = None
        
        if drugs[drug]["all"]:
            # COBRAモデル内の遺伝子の有無を判定
            drugs[drug]["cobra"], another = checkinCOBRA(drugs[drug]["all"], model_c)
            if drugs[drug]["cobra"]:
                print " >>> %s has metabolic target" % drug
                for gene in drugs[drug]["cobra"]:
                    if not gene in cobra_target:
                       cobra_target[gene] += drugs[drug]["a_ex"]

            # riboモデル内の判定
            # riboはシングルのみ(複数の場合もデフォルトで固定)
            for gene in another:
                if gene in ribo_bnum: # ribo遺伝子が含まれている場合
                    ribo_data["flag"] = True
                    ribo_data["a_ex"] += drugs[drug]["a_ex"] * 0.369 # riboのDose調節
                    if dataset_r.get(drug): # 薬剤のriboデータがある場合
                        ribo_data["dataset"] = dataset_r[drug]["data"]
                        print " >>> %s has ribosome target" % drug

    # run phase (フィードバックはまだ考えない)
    result_r_list = []
    ribo_res = 1 # riboのλ の返り値
    # 1回目から始めると、wt側に入らなくなる
    count = 1 # ループの回数
    before = None # ループ中の前回の結果

    while 1: # while文にする
        # COBRA
        model_c = wt_model.copy() # modelの初期化
        flux = wt_flux.copy() # fluxの初期化

        if count == 0 or not cobra_target:
            Lambda_0 = model_c.solution.f # growth rate of wt
            count += 1
        else: # COBRAをターゲットにしている薬剤の分繰り返す
            cobra_target = {gene: 1.0 / (cobra_target[gene] + 1.0) 
                            for gene in cobra_target.keys()}
            result_c = multi_knockdown(model_c, cobra_target,
                                       wt_flux=flux)
            Lambda_0 = result_c[1].values()[0] # growth rate of knockdown

        # 終了判定
        if before and abs(before - Lambda_0) <= 0.001 * abs(before):
            break

        ribo_data["dataset"]["Lambda_0"] = Lambda_0

        # ribo
        # singleのみ
        if ribo_data["flag"]:
            result_ribo = run(ribo_data["a_ex"], ribo_data["dataset"], step=step_int)
            ribo_data["dataset"] = result_ribo["dataset"]
            ribo_res = result_ribo["result"]
            ribo_data["a_ex"] = result_ribo["a_ex"]
   
            result_r_list.append(result_ribo["result"]) # allDataの方がいいかも?
        else:
            ribo_res = Lambda_0 # riboに入らなかった場合は

        before = ribo_res
        
        # 1回だけ
        break
        
    # result_r_listにriboモデルの計算結果たちが並んでいる
    # return result_c, result_r_list
    # return Lambda_0, result_r_list
    return ribo_res, (Lambda_0, result_r_list)

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pylab as plt

    result = []
    Lambda_0 = 0
    
    for dose in np.linspace(0, 1., 201): 
        growth, (cobra, ribo) = crModel({"DB01332": dose, "DB01082": dose})
        result.append([dose, growth])
        Lambda_0 = cobra
        print Lambda_0

    result = np.array(result)
    plt.plot(result.T[0], result.T[1], "og")
    plt.xlabel(r"Dose")
    plt.ylabel(r"Growth Rate $\lambda/\lambda_{0}$")

    plt.savefig("result/cobra.png", dpi=200)
    
    with open("result/cobra.csv", "w") as wf:
        for i in result:
            wf.write("%e, %e\n" % (i[0], i[1]))



