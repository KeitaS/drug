#coding:utf-8

import numpy as np
import matplotlib.pylab as plt
from ecell4 import *
util.decorator.SEAMLESS_RATELAW_SUPPORT = True


def createModel(r_max=65.8, r_min=19.3, K_D=1.0, K_t=6.1*10**-2, K_on=3.0, Lambda_0=1.35, Lambda_0_a=0.31, IC50=0.41, IC50_a=0.189):
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

    薬剤によって変更する必要がある値:
    Lambda_0_a, IC50, IC50_a

    培地によって変更する値(COBRAの結果):
    Lambda_0

    """

    Delta_r = r_max - r_min # µM
    K_off = K_on * K_D

    P_in = Delta_r * Lambda_0_a / 2.0 / IC50_a
    P_out = (Lambda_0_a / 2) ** 2.0 / K_t / K_D

    with reaction_rules():
        Lambda = (r_u - r_min) * K_t

        ~a_ex > a | P_in * a_ex
        a > ~a_ex | P_out * a
        a + r_u > r_b | K_on * a * (r_u - r_min)
        r_b > a + r_u | K_off * r_b
        ~r_u > r_u | Lambda * (r_max - Lambda * Delta_r * (1 / Lambda_0 - 1/K_t / Delta_r))
        a > ~a | a * Lambda
        r_u > ~r_u | r_u * Lambda
        r_b > ~r_b | r_b * Lambda

    return get_model()

def run_test2(a_ex, step=10., legend=[], inpData={}, y0={"a": .0, "r_u": 30.0, "r_b": .0}):
    dataset = {"Lambda_0": 1.35, "Lambda_0_a": 0.31, "IC50": 0.41, "IC50_a": 0.189, "K_t": 6.1 * 10 ** -2, "r_min": 19.3}

    dataset.update(inpData)
    model = createModel(**dataset)
    y0["a_ex"] = a_ex

    if not legend:
        legend = y0.keys()

    runsim = run_simulation(step, solver="ode", y0=y0,
                            return_type="observer", model=model,
                            species_list=legend)
    data = runsim.data()
    return data, legend


def makeGraph(data, savename, legend=[], title="", xlabel="", ylabel=""):
    for i in range(len(data[0]) - 1):
        if legend[i]:
            plt.plot(data.T[0], data.T[i+1], label=legend[i])
        else:
            plt.plot(data.T[0], data.T[i+1])
    if title: # titleがあったらつける
        plt.title(title)

    if xlabel: # xlabelがあったらつける
        plt.xlabel(xlabel)

    if ylabel: # ylabelがあったらつける
        plt.ylabel(ylabel)

    plt.legend(loc="upper right")
    # plt.show()
    plt.savefig("result/%s" % (savename), dpi=200)
    plt.close()


def run_test(dataset={}, y0={"a": .0, "r_u": 30.0, "r_b": .0}, step=[0, 1], stepInt=200):
    """
    riboモデルをRunするモジュール
    現段階では、Lambda

    y0: 各変数の初期値{"a": .0, "r_u": 30.0, "r_b": .0}, a_exは振る。
    dataset: riboモデルに入力必須のもの。Lambda_0, Lambda_0_a, IC50, IC50_a。
    step: a_exの振り幅。default = [0, 1]
    # streptmycin, canamycin: 1, tetracycline: 2, chloramphenicol: 20
    """

    default_data = {"Lambda_0": 1.35, "Lambda_0_a": 0.31, "IC50": 0.41, "IC50_a": 0.189, "K_t": 6.1 * 10 ** -2, "r_min": 19.3}

    for variable in default_data.keys():
        if not dataset.get(variable):
            dataset[variable] = default_data[variable]

    model = createModel(**dataset)

    result = []
    for a_ex in np.linspace(step[0], step[1], 201):
        y0["a_ex"] = a_ex
        runsim = run_simulation((0, stepInt), solver = "ode", y0 = y0,
                                return_type = "observer", model = model,
                                species_list = ["r_u"])
        r_u = runsim.data()[-1][1]
        Lambda = (r_u - dataset["r_min"]) * dataset["K_t"]
        result.append((a_ex, Lambda / dataset["Lambda_0"]))

    return np.array(result)



def run(a_ex, dataset={}, y0={"a": .0, "r_u": 30.0, "r_b": .0}, step=10):
    """
    12/16
    riboモデルをRunするモジュール
    a_exも入力するようになっている。
    ### sim_time: 必要なら入れる ###
    dataset: riboモデルに入力必須のもの。Lambda_0, Lambda_0_a, IC50, IC50_a。
    y0: 各変数の初期値{"a": .0, "r_u": 30.0, "r_b": .0}, a_exは振る。
    step: a_exの振り幅。default = 100
    # streptmycin, canamycin: 1, tetracycline: 2, chloramphenicol: 20
    """

    default_data = {"Lambda_0": 1.35, "Lambda_0_a": 0.31, "IC50": 0.41, "IC50_a": 0.189, "K_t": 6.1 * 10 ** -2, "r_min": 19.3}

    #for variable in default_data.keys():
    #    if not dataset.get(variable):
    #        dataset[variable] = default_data[variable]
    default_data.update(dataset)

    model = createModel(**default_data)

    result = {}

    y0["a_ex"] = a_ex
    runsim = run_simulation(step, solver="ode", y0=y0,
                            return_type="observer", model=model,
                            species_list=["a", "a_ex", "r_u", "r_b"])

    data = runsim.data()[-1]
    r_u = data[3]
    Lambda = (r_u - default_data["r_min"]) * default_data["K_t"]
    result["a_ex"] = data[2]
    result["dataset"] = {"a": data[1], "r_u": r_u, "r_b": data[4]}
    # result["result"] = (sim_time + data[0], Lambda / dataset["Lambda_0"])
    result["growth"] = Lambda # Lambdaの結果を返すように
    return result


if __name__ == "__main__":
    r_min = 19.3
    K_t = 6.1 * 10 ** -2

    ## drug data
    drug = ["streptmycin", "kanamycin", "tetracycline", "chloramphenicol"]
    Lambda_0 =  [1.35, 0.85, 0.40] # 1/h (1.35, 0.85, 0.40)
    Lambda_0_a = {"streptmycin": 0.31, "kanamycin": 0.169, "tetracycline": 5.24, "chloramphenicol": 1.83} # 1/h
    IC50 = {"streptmycin": [0.41, 0.28, 0.196], "kanamycin": [0.246, 0.096, 0.065], "tetracycline": [0.5, 0.6, 1.45], "chloramphenicol": [2.85, 2.65, 5.7]} # µg/ml
    IC50_a = {"streptmycin": 0.189, "kanamycin": 0.05, "tetracycline": 0.229, "chloramphenicol": 2.49} # µg/ml
    A_ex = {"streptmycin": 0.6, "kanamycin": 0.5, "tetracycline":2, "chloramphenicol": 20}

    xlabel = "Extracellular antibiotic concentration $a_{ex}$ ($\mu$M)"
    ylabel = "Normalized Growth Rate $\lambda/\lambda_{0}$"

    for name in drug:
        print "%s simulation >>>" % (name)
        data = []
        for i in range(3):
            print "medium %d >>" % (i)
            dataset = {"Lambda_0": Lambda_0[i],
                       "Lambda_0_a": Lambda_0_a[name],
                       "IC50": IC50[name][i],
                       "IC50_a":IC50_a[name]}
            legend = ["r_u"]
            count = 0
            for j in np.linspace(0, A_ex[name], 51):
                print count
                result, legend = run_test2(j, inpData=dataset, legend=legend)

                if i == 0:
                    result = [j, (result[-1][1] - r_min) * K_t / Lambda_0[i]]
                    data.append(result)
                else:
                    data[count].append((result[-1][1] - r_min) * K_t / Lambda_0[i])
                count += 1

        savename = "20160607/4/%s.png" % (name)
        legend = ["$Gly$", "$Gly_{CAA}$", "$Gly_{RDM}$"]
        makeGraph(np.array(data), savename, legend, name, xlabel, ylabel)
