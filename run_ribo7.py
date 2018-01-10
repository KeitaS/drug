# coding: utf-8
from modules import *

if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    
    # combinatorial simulation
    ## merge DataFiles
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # csvdir = "results/ribo7/csv/sim100"
    # for drugName in drugNameList:
    #     dirName = "{}/{}".format(csvdir, "_".join(drugName))
    #     fileNameList = ["{}/{}_{}.csv".format(dirName, "_".join(drugName), num) for num in range(101)]
    #     df = mergeResults(fileNameList)
    #     df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(drugName)), index=False)

    ## SameDrug combination
    drugNameList = [[name, name] for name in drugNames]
    dataLists = ["results/ribo7/csv/sim100/{}_merge.csv".format("_".join(i)) for i in drugNameList]
    doubleSaveName = "results/ribo7/images/sameDrug_sim100.png"
    createHeatmap(drugNameList, dataLists, [2, 2], doubleSaveName)


    ## differentDrug combination
    drugNameList = list(itr.combinations(drugNames, 2))
    dataLists = ["results/ribo7/csv/sim100/{}_merge.csv".format("_".join(i)) for i in drugNameList]
    doubleSaveName = "results/ribo7/images/diffDrug_sim100.png"
    createHeatmap(drugNameList, dataLists, [3, 2], doubleSaveName)

    # combinatorial simulation(virtual drug)
    ## merge DataFiles
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["r30", "r30"], ["r30", "r50"]]
    # csvdir = "results/ribo7/csv/sim100_v"
    # for drugName in drugNameList:
    #     for target in targetList:
    #         dirName = "{}/{}".format(csvdir, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
    #         fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
    #         df = mergeResults(fileNameList)
    #         df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))])), index=False)
            
    ## create Image
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["r30", "r30"], ["r30", "r50"]]
    nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for target in targetList for drugName in drugNameList]
    dataLists = ["results/ribo7/csv/sim100_v/{}_merge.csv".format("_".join(name)) for name in nameList]
    saveName = "results/ribo7/images/virtualDrug_sim100.png"
    createHeatmap(nameList, dataLists, [3, 2], saveName)

    # oldeval simulation
    ## merge DataFiles
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # csvdir = "results/ribo7/csv/old100"
    # for drugName in drugNameList:
    #     fileNameList = ["{}/{}/{}.csv".format(csvdir, "_".join(drugName), num) for num in range(101)]
    #     df = mergeResults(fileNameList)
    #     df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(drugName)), index=False)
    
    ## SameDrug Combination
    drugNameList = [[name, name] for name in drugNames]
    csvdir = "results/ribo7/csv/old100"
    dataList = ["{}/{}_merge.csv".format(csvdir, "_".join(i)) for i in drugNameList]
    saveName = "results/ribo7/images/sameDrug_old100.png"
    createHeatmap(drugNameList, dataList, [2, 2], saveName, ["a2", "a1"], "epsilon")

    ## DiffDrug Combination
    drugNameList = list(itr.combinations(drugNames, 2))
    csvdir = "results/ribo7/csv/old100"
    dataList = ["{}/{}_merge.csv".format(csvdir, "_".join(i)) for i in drugNameList]
    saveName = "results/ribo7/images/diffDrug_old100.png"
    createHeatmap(drugNameList, dataList, [3, 2], saveName, ["a2", "a1"], "epsilon")

    # oldeval simulation (virtual drug)
    ## merge DataFiles
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["r30", "r30"], ["r30", "r50"]]
    # csvdir = "results/ribo7/csv/old100_v"
    # nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for target in targetList for drugName in drugNameList]
    # for name in nameList:
    #     dirName = "{}/{}".format(csvdir, "_".join(name))
    #     fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
    #     df = mergeResults(fileNameList)
    #     df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(name)), index=False)

    ## create Image
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["r30", "r30"], ["r30", "r50"]]
    csvdir = "results/ribo7/csv/old100_v"
    nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for target in targetList for drugName in drugNameList]
    dataList = ["{}/{}_merge.csv".format(csvdir, "_".join(name)) for name in nameList]
    saveName = "results/ribo7/images/virtualDrug_old100.png"
    createHeatmap(nameList, dataList, [3, 2], saveName, ["a2", "a1"], "epsilon")

    # neweval simulation
    ## merge DataFiles
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # csvdir = "results/ribo7/csv/new100"
    # for drugName in drugNameList:
    #     fileNameList = ["{}/{}/{}.csv".format(csvdir, "_".join(drugName), num) for num in range(5)]
    #     df = mergeResults(fileNameList)
    #     df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(drugName)), index=False)

    ## SameDrug Combination
    drugNameList = [[name, name] for name in drugNames]
    csvdir = "results/ribo7/csv/new100"
    dataList = ["{}/{}_merge.csv".format(csvdir, "_".join(i)) for i in drugNameList]
    saveName = "results/ribo7/images/sameDrug_new100.png"
    createHeatmap(drugNameList, dataList, [2, 2], saveName, ["Slope", "MidPoint"], "LinerType")

    ## DiffDrug Combination
    drugNameList = list(itr.combinations(drugNames, 2))
    csvdir = "results/ribo7/csv/new100"
    dataList = ["{}/{}_merge.csv".format(csvdir, "_".join(i)) for i in drugNameList]
    saveName = "results/ribo7/images/diffDrug_new100.png"
    createHeatmap(drugNameList, dataList, [3, 2], saveName, ["Slope", "MidPoint"], "LinerType")

    # neweval simulation (virtual drug)
    ## merge DataFiles
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["r30", "r30"], ["r30", "r50"]]
    # csvdir = "results/ribo7/csv/new100_v"
    # nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for drugName in drugNameList for target in targetList]
    # for name in nameList:
    #     dirName = "{}/{}".format(csvdir, "_".join(name))
    #     fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(5)]
    #     df = mergeResults(fileNameList)
    #     df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(name)), index=False)

    ## createImage
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["r30", "r30"], ["r30", "r50"]]
    csvdir = "results/ribo7/csv/new100_v"
    nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for target in targetList for drugName in drugNameList]
    dataList = ["{}/{}_merge.csv".format(csvdir, "_".join(i)) for i in nameList]
    saveName = "results/ribo7/images/virtualDrug_new100.png"
    createHeatmap(nameList, dataList, [3, 2], saveName, ["Slope", "MidPoint"], "LinerType")

   

