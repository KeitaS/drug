#coding:utf-8

"""
DrugBankのデータを遺伝子番号に変更し、オブジェクトにして返すモジュール。
key: 薬剤名
value: b-numberが入ったリスト
"""

import re, os, collections, copy


def drug2upid(fname="model/approved_target_ids_all.csv"):
    """
    DrubBankから取ってきたデータを抽出するモジュール
    result: key = Drug name
            value = list of UniPlot Id
    """
    with open(fname, "r") as fopen:
        result = collections.defaultdict(list)
        for index, line in enumerate(fopen.readlines()):
            line = re.split(",", line.strip())
            if index == 0:
                target = ["UniProt ID", "Species", "Drug IDs"]
                for num, title in enumerate(line):
                    if title in target:
                        target.append(num)
                del target[0:3]
            
            else:
                # drugData = [upid, species, drug ids]
                drugData = [line[target[x]] for x in range(len(target))]
                if re.search("Escherichia coli", drugData[1]):
                    drugList = re.split("; ", drugData[2])
                    for drug in drugList:
                        result[drug].append(drugData[0])

    return result


def upid2bnum(data):
    """
    drug2upidで作成したデータを用いて、UniPlotIdをb-numberに変更するモジュール
    """
    ass_dict = association()
    for key, value in data.items():
        for index, upid in enumerate(value):
            if ass_dict.get(upid):
                data[key][index] = ass_dict.get(upid)

    return data
            

def association(fname="model/ecoli.txt"):
    """
    UniPlotIdとb-numberの対応表をdict型にするモジュール
    """
    ass_dict = {}
    with open(fname, "r") as fopen:
        for line in fopen.readlines():
            if re.match("^b[0-9]*", line):
                line = re.split("[;\s\t]*", line.strip())
                ass_dict[line[3]] = line[0]

    return ass_dict


    
def drug2bnum(fname="model/approved_target_ids_all.csv"):
    data = drug2upid(fname)
    result = upid2bnum(data)
    
    return result


if __name__ == "__main__":
    import sys
    data = drug2bnum("model/approved_target_ids_all.csv")
    print data.keys()
    print len(data.keys())
