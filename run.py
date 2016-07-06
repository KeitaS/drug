#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True

@reaction_rules
def r30_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a > ~a | a * Lambda
    a + r30_u > r30_b | K_on * a * (r_u - r_min)
    r30_b > a + r30_u | K_off * r30_b

@reaction_rules
def r50_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a > ~a | a * Lambda
    a + r50_u > r50_b | K_on * a * (r_u - r_min)
    r50_b > a + r50_u | K_off * r50_b

@reaction_rules
def ribo_binding_reaction(a_ex, a, P_in, P_out, K_on, K_off, Lambda):
    ~a_ex > a | P_in * a_ex
    a > ~a_ex | P_out * a
    a > ~a | a * Lambda
    a + r_u > r_b | K_on * a * (r_u - r_min)
    r_b > a + r_u | K_off * r_b
    r_b > a + r30_u + r50_u | Kd * r_u # r_bの解離


def createModel(drugs=[], r_max=65.8, r_min=19.3, K_D=1.0, K_t=6.1*10**-2, K_on=3.0, Lambda_0=1.35, Kd=100., p=1., target=[]):
    """
    リボソームモデルを構成するモジュール
    r_max: µM
    r_min: µM
    K_D: none
    K_t: 1/µM/h
    K_on: ???
    K_off: ???
    Delta_r: µM
    Lambda_0: default medium = Gly_RDM
    Lambda_0_a: default drug = streptmycin
    IC50: default medium = Gly_RDM, default drug = streptmycin
    IC50_a: default drug = streptmycin
    Ka: subunitの結合
    Kd: subunitの解離

    薬剤によって変更する必要がある値:
    Lambda_0_a, IC50, IC50_a

    培地によって変更する値(COBRAの結果):
    Lambda_0

    Consideration:
    薬剤固有の値はIC50_a, IC50, Lambda_0_a
    複数入ってきた時の対応はどうするか
    """

    Delta_r = r_max - r_min # µM
    K_off = K_on * K_D # riboと薬剤との結合
    r_u_0 = Lambda_0 / K_t + r_min
    Ka = (Lambda_0 + Kd) / (p * p * r_u_0)

    # P_in, P_outをdrugのデータに入れる
    for drug in drugs: # drug dataに基づく情報を用いた
        drug["P_in"] = Delta_r * drug["Lambda_0_a"] / 2.0 / drug["IC50_a"] # 薬剤の流入
        drug["P_out"] = (drug["Lambda_0_a"] / 2) ** 2.0 / K_t / K_D # 薬剤の流出


    with reaction_rules():
        ### expression
        Lambda = (r_u - r_min) * K_t
        SUP = (Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))) * (1 + p) # subunit production expression

        ### reaction
        ## drug
        if len(drugs) > 0:
            if drugs[0]["type"] == "30s":
                # print "drug1 targets 30s ribosomal subunit >>"
                r30_binding_reaction(a1_ex, a1, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
            elif drugs[0]["type"] == "50s":
                # print "drug1 targets 50s ribosomal subunit >>"
                r50_binding_reaction(a1_ex, a1, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
            elif drugs[0]["type"] == "ribo":
                # print "drug1 targets ribosome >>"
                ribo_binding_reaction(a1_ex, a1, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)

        if len(drugs) > 1:
            if drugs[1]["type"] == "30s":
                # print "drug2 targets 30s ribosomal subunit >>"
                r30_binding_reaction(a2_ex, a2, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
            elif drugs[1]["type"] == "50s":
                # print "drug2 targets 50s ribosomal subunit >>"
                r50_binding_reaction(a2_ex, a2, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)
            elif drugs[1]["type"] == "ribo":
                # print "drug2 targets ribosome >>"
                ribo_binding_reaction(a2_ex, a2, drugs[0]["P_in"], drugs[0]["P_out"], K_on, K_off, Lambda)

        ## ribo and subunit
        # production
        ~r30_u > r30_u | SUP
        ~r50_u > r50_u | SUP

        #bonding
        r30_u + r50_u > r_u | Ka * r30_u * r50_u

        # dissociation
        r_u > r30_u + r50_u | Kd * r_u
        r_u > ~r_u | r_u * Lambda
        r30_u > ~r30_u | r30_u * Lambda
        r50_u > ~r50_u | r50_u * Lambda

        # diffusion
        r30_b > ~r30_b | r30_b * Lambda
        r50_b > ~r50_b | r50_b * Lambda
        r_b > ~r_b | r_b * Lambda

    return get_model()


def checkRiboDrug(drugNames=[]):
    """
    薬剤の標的を返す関数
    """
    targets = []
    for drugName in drugNames:
        if drugName == "30s":
            targets.append("30s")
        elif drugName == "50s":
            targets.append("50s")
        elif drugName == "ribo":
            targets.append("ribo")

    return targets


def run(drugs=[], step=10., legend=[], inpData={}, y0={"r30_u": 30., "r50_u": 30., "r_u": 30., "r_b": .0}):
    """
    ribosomeモデル実行関数

    drugs: [{"name":, "type":, "dose":, }]
    step:
    legend:
    inpData:
    y0:

    ToDo:
    今後、drugの種類を判定する関数をつくり、この関数に導入する。
    """

    dataset = {"Lambda_0": 1.35, "K_t": 6.1 * 10 ** -2, "r_min": 19.3}

    dataset.update(inpData) # inpDataの内容をdatasetに入れる
    if drugs:
        dataset["drugs"] = drugs
        # y0に薬剤のdoseを入れる
        for index, drug in enumerate(drugs):
            y0["a%d_ex" % (index + 1)] = drug["dose"]

    model = createModel(**dataset) # モデルを作成

    # legend
    if not legend:
        legend = y0.keys()

    runsim = run_simulation(step, solver="ode", y0=y0,
                            return_type="observer", model=model,
                            species_list=legend)
    data = runsim.data()
    return data, legend


def makeGraph(data, savename, legend=[], title="", xlabel="", ylabel=""):
    """
    グラフ作成用モジュール
    """
    for i in range(len(data[0]) - 1):
        if len(legend) > i:
            plt.plot(data.T[0], data.T[i+1], label=legend[i])
        else:
            plt.plot(data.T[0], data.T[i+1])

    if title: # titleがあったらつける
        plt.title(title)

    if xlabel: # xlabelがあったらつける
        plt.xlabel(xlabel)

    if ylabel: # ylabelがあったらつける
        plt.ylabel(ylabel)

    if legend:
        plt.legend(loc="upper right")
    # plt.show()
    plt.savefig(savename, dpi=200)
    plt.close()


if __name__ == "__main__":
    r_min = 19.3
    K_t = 6.1 * 10 ** -2
    Kd = 5.
    p = 0.1

    ## drug data
    drug = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    Lambda_0 =  [1.35, 0.85, 0.40] # 1/h (1.35, 0.85, 0.40)
    Lambda_0_a = {"Streptmycin": 0.31, "Kanamycin": 0.169, "Tetracycline": 5.24, "Chloramphenicol": 1.83} # 1/h
    IC50 = {"Streptmycin": [0.41, 0.28, 0.196], "Kanamycin": [0.246, 0.096, 0.065], "Tetracycline": [0.5, 0.6, 1.45], "Chloramphenicol": [2.85, 2.65, 5.7]} # µg/ml
    IC50_a = {"Streptmycin": 0.189, "Kanamycin": 0.05, "Tetracycline": 0.229, "Chloramphenicol": 2.49} # µg/ml
    a_ex = {"Streptmycin": 0.6, "Kanamycin": 0.5, "Tetracycline":2, "Chloramphenicol": 20}

    ## Kanamycinで比較
    name = "Kanamycin"
    dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}
    legend = ["r_u"]

    # 保存用ディレクトリの作成
    import os
    if not os.path.exists("images/result"):
        os.makedirs("images/result")
    del(os)

    # 以前のモデルと今回のモデルの比較
    ## ribo3
    """
    drugs = [{"name": name, "type": "ribo", "dose": .0, "Lambda_0_a": Lambda_0_a[name], "IC50": IC50[name][0], "IC50_a": IC50_a[name]}]
    drug_types = ["30s", "50s", "ribo"]
    single_result = [] # 単剤用の結果の格納場所
    for num in range(len(drug_types)):
        drug_type = drug_types[num]
        drugs[0]["type"] = drug_type

        # doseを振り、モデルを実行
        for count, dose in enumerate(np.linspace(0, a_ex[name], 51)):
            drugs[0]["dose"] = dose
            result, legend = run(drugs, inpData=dataset, legend=legend)
            result = (result[-1][1] - r_min) * K_t / Lambda_0[0] # 結果をgrowthに書き換え
            if num == 0:
                single_result.append([dose, result])
            else:
                single_result[count].append(result)


    ## ribo.py
    from ribo import run as run1
    dataset = {"Lambda_0": Lambda_0[0], "Lambda_0_a": Lambda_0_a[name], "IC50": IC50[name][0], "IC50_a":IC50_a[name]}
    for index, dose in enumerate(np.linspace(0, a_ex[name], 51)):
        result = run1(dose, dataset=dataset)
        single_result[index].append(result["growth"] / Lambda_0[0])

    print single_result

    ## make Graph
    legend = ["30s", "50s", "ribo", "before"]
    savename = "images/result/ribo_ribo3.png"
    xlabel = "Extracellular Antibiotic Concentration $a_{ex}$ ($\mu$M)"
    ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"
    title = "Single Drug Reaction"
    makeGraph(np.array(single_result), savename, legend, title, xlabel, ylabel)
    """


    # nature2006のモデルとの比較
    import matplotlib.pylab as plt
    dataset = {"Kd": Kd, "Lambda_0": Lambda_0[0], "p": p}

    ## Chloramphenicol 分子量: 323.132 g/mol, Dose: 1 µg/ml
    Chloramphenicol = 1. * 10 ** 3 / 323.132 # µM
    ## Tetracycline 分子量: 444.435 g/mol, Dose: 2 µg/ml
    Tetracycline = 2. * 10 ** 3 / 444.435 #µM
    drug_data = {"Streptmycin": {"dose": 5., "type": "30s"},
                 "Chloramphenicol": {"dose": Chloramphenicol, "type": "50s"},
                 "Tetracycline": {"dose": Tetracycline, "type": "30s"}
                }
    print drug_data

    drug_comb = [["Chloramphenicol", "Tetracycline"],
                 ["Chloramphenicol", "Streptmycin"],
                 ["Tetracycline", "Streptmycin"]
                ]
    legend = ["r_u"]

    for index, drug in enumerate(drug_comb):
        # φ
        second_result = [1.0]
        name1 = drug[0]
        name2 = drug[1]

        # 1つめ
        drugs = [{"name": name1, "type": drug_data[name1]["type"], "dose": drug_data[name1]["dose"], "Lambda_0_a": Lambda_0_a[name1], "IC50": IC50[name1][0], "IC50_a": IC50_a[name1]}]
        result, legend = run(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
        second_result.append(result)

        # 2つめ
        drugs = [{"name": name2, "type": drug_data[name2]["type"], "dose": drug_data[name2]["dose"], "Lambda_0_a": Lambda_0_a[name2], "IC50": IC50[name2][0], "IC50_a": IC50_a[name2]}]
        result, legend = run(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
        second_result.append(result)

        # 組み合わせ
        drugs = [{"name": name1, "type": drug_data[name1]["type"], "dose": drug_data[name1]["dose"], "Lambda_0_a": Lambda_0_a[name1], "IC50": IC50[name1][0], "IC50_a": IC50_a[name1]},
                 {"name": name2, "type": drug_data[name2]["type"], "dose": drug_data[name2]["dose"], "Lambda_0_a": Lambda_0_a[name2], "IC50": IC50[name2][0], "IC50_a": IC50_a[name2]}
                ]
        result, legend = run(drugs, inpData=dataset, legend=legend)
        result = (result[-1][1] - r_min) * K_t / Lambda_0[0]
        second_result.append(result)

        # グラフ化
        x = range(4)
        plt.subplot(330 + index + 4)
        plt.bar(x, second_result, align="center")
        plt.xticks(x, ["$\phi$", name1[:3], name2[:3], "%s+\n%s" % (name1[:3], name2[:3])])
        plt.xlabel("Drug Name")
        plt.ylabel("Normalized Growth Rate $\lambda/\lambda_{0}$")
        plt.ylim(0, 1.1)
        plt.title("%s & %s" % (name1[:3], name2[:3]))
    plt.tight_layout()
    plt.savefig("images/result/double_bar.png", bbox_inches="tight", pad_inches=0.3, dpi=200)
    plt.close()