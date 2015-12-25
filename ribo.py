#coding:utf-8

import numpy as np
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



def run(a_ex, dataset={}, y0={"a": .0, "r_u": 30.0, "r_b": .0}, step=100):
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

    for variable in default_data.keys():
        if not dataset.get(variable):
            dataset[variable] = default_data[variable]

    model = createModel(**dataset)
    
    result = {}

    y0["a_ex"] = a_ex
    runsim = run_simulation((0, step), solver="ode", y0=y0, 
                            return_type="observer", model=model,
                            species_list=["a", "a_ex", "r_u", "r_b"])
    
    data = runsim.data()[-1]
    r_u = data[3]
    Lambda = (r_u - dataset["r_min"]) * dataset["K_t"]
    result["a_ex"] = data[2]
    result["dataset"] = {"a": data[1], "r_u": r_u, "r_b": data[4]}
    # result["result"] = (sim_time + data[0], Lambda / dataset["Lambda_0"])
    result["result"] = Lambda # Lambdaの結果を返すように
    return result


if __name__ == "__main__":
    dataset = {"Lambda_0": 0.982371812727}
    result = run(5.0, dataset=dataset) 
    print result

