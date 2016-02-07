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
    $B$"$k$b$N$H$J$$$b$N$I$A$i$b(Breturn$B$9$k(B
    cobra$B$KF~$C$F$$$k0dEA;R$N8!:wMQ%b%8%e!<%k(B

    genes: $B8!:w$7$?$$0dEA;R$N%j%9%H(B
    model: COBRA$B$G:n@.$7$?%b%G%k(B
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
    drug id$B$+$i=i4|CMEy$N%G!<%?$rCj=P(B
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
    cobra$B$H(Bribo$B$NJ#9g%b%G%k(B
    COBRA$B$NH=JL$rL\;X$9(B
    drugs$B$G(Ba_ex$B$bEO$;$k$h$&$K(B

    drugs: {DrugBank ID: Dose}
    step_int: riboModel step
    model_c: COBRA model
    model_r: riboModel data
    drug_data: drug data was downloaded from DrugBank. 

    ribo$B$N(BDose$BD4@0$K$O(Bstreptmycin

    ToDo:
        ribo$B$N(BDose$BD4@0$rA4Lt:^$GE,1~(B
        ribo$B$N(Binput$B$r>\:Y2=(B
        ribo$B$+$i(Bcobra$B$N%U%#!<%I%P%C%/(B
    """

    # load phase
    model_c = read_sbml_model(model_c)
    model_c.optimize(solver="glpk")
    wt_flux = model_c.solution.x_dict.copy()
    wt_model = model_c.copy()
    drugs = {drug:{"a_ex": drugs[drug]} for drug in drugs.keys()} # drug$BL>$r(Bkey$B$K$7$?(Bdict$B$K$J$k(B
    ribo_bnum = ["b3230", "b3319", "b3342", "b3296", "b3298", "P61177", "b3307", "P60725", "b3316", "b3313"] # ribo$B$r%?!<%2%C%H$H$7$F$$$k(Bbnumber, $B$3$3$rA}$d$7$F$$$/(B.

    # check phase
    drug_data = drug2bnum(drug_data)
    dataset_r = makeRiboData(model_r) # {drug_id:{"data":{"Lambda_0_a":, "IC50":, "IC50_a":,}, drug name, b-number}}
    
    cobra_target = defaultdict(lambda: 0.0) # cobra$B$r%?!<%2%C%H$K$9$k0dEA;R$H(Bfold_change$B$N(Bdict
    ribo_target = [] # ribo$B$r%?!<%2%C%H$K$9$kLt:^$N%j%9%H(B

    # $B%j%\%=!<%`$OC10lLt:^$N$_$J$N$G!"%?!<%2%C%H$K$J$C$F$$$?$i(Bfrag$B$r(BTrue$B$K$9$k(B
    ribo_data = {"flag": False, "a_ex": 0, "dataset": {}}
    
    for drug in drugs.keys():
        if drug_data.get(drug):
            drugs[drug]["all"] = drug_data[drug]
        else:
            drugs[drug]["all"] = None

        drugs[drug]["cobra"] = None
        drugs[drug]["ribo"] = None
        
        if drugs[drug]["all"]:
            # COBRA$B%b%G%kFb$N0dEA;R$NM-L5$rH=Dj(B
            drugs[drug]["cobra"], another = checkinCOBRA(drugs[drug]["all"], model_c)
            if drugs[drug]["cobra"]:
                print " >>> %s has metabolic target" % drug
                for gene in drugs[drug]["cobra"]:
                    if not gene in cobra_target:
                       cobra_target[gene] += drugs[drug]["a_ex"]

            # ribo$B%b%G%kFb$NH=Dj(B
            # ribo$B$O%7%s%0%k$N$_(B($BJ#?t$N>l9g$b%G%U%)%k%H$G8GDj(B)
            for gene in another:
                if gene in ribo_bnum: # ribo$B0dEA;R$,4^$^$l$F$$$k>l9g(B
                    ribo_data["flag"] = True
                    ribo_data["a_ex"] += drugs[drug]["a_ex"] * 0.369 # ribo$B$N(BDose$BD4@a(B
                    if dataset_r.get(drug): # $BLt:^$N(Bribo$B%G!<%?$,$"$k>l9g(B
                        ribo_data["dataset"] = dataset_r[drug]["data"]
                        print " >>> %s has ribosome target" % drug

    # run phase ($B%U%#!<%I%P%C%/$O$^$@9M$($J$$(B)
    result_r_list = []
    ribo_res = 1 # ribo$B$N&K(B $B$NJV$jCM(B
    # 1$B2sL\$+$i;O$a$k$H!"(Bwt$BB&$KF~$i$J$/$J$k(B
    count = 1 # $B%k!<%W$N2s?t(B
    before = None # $B%k!<%WCf$NA02s$N7k2L(B

    while 1: # while$BJ8$K$9$k(B
        # COBRA
        model_c = wt_model.copy() # model$B$N=i4|2=(B
        flux = wt_flux.copy() # flux$B$N=i4|2=(B

        if count == 0 or not cobra_target:
            Lambda_0 = model_c.solution.f # growth rate of wt
            count += 1
        else: # COBRA$B$r%?!<%2%C%H$K$7$F$$$kLt:^$NJ,7+$jJV$9(B
            cobra_target = {gene: 1.0 / (cobra_target[gene] + 1.0) 
                            for gene in cobra_target.keys()}
            result_c = multi_knockdown(model_c, cobra_target,
                                       wt_flux=flux)
            Lambda_0 = result_c[1].values()[0] # growth rate of knockdown

        # $B=*N;H=Dj(B
        if before and abs(before - Lambda_0) <= 0.001 * abs(before):
            break

        ribo_data["dataset"]["Lambda_0"] = Lambda_0

        # ribo
        # single$B$N$_(B
        if ribo_data["flag"]:
            result_ribo = run(ribo_data["a_ex"], ribo_data["dataset"], step=step_int)
            ribo_data["dataset"] = result_ribo["dataset"]
            ribo_res = result_ribo["result"]
            ribo_data["a_ex"] = result_ribo["a_ex"]
   
            result_r_list.append(result_ribo["result"]) # allData$B$NJ}$,$$$$$+$b(B?
        else:
            ribo_res = Lambda_0 # ribo$B$KF~$i$J$+$C$?>l9g$O(B

        before = ribo_res
        
        # 1$B2s$@$1(B
        break
        
    # result_r_list$B$K(Bribo$B%b%G%k$N7W;;7k2L$?$A$,JB$s$G$$$k(B
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



