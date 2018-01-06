#coding: utf-8
from run_modules import *

if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    csvdir = "results/ribo5/csv/sim100"
    
    # combinatorial simulation
    ## merge DataFiles
    # for modif in range(1, 3):
    #     for drugName in drugNameList:
    #         dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(drugName))
    #         fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
    #         df = mergeResults(fileNameList)
    #         df.to_csv("{}/modif{}/{}_merge.csv".format(csvdir, modif, "_".join(drugName)), index=False)

    ## SameDrug combination
    # drugNameList = [[name, name] for name in drugNames]
    # for modif in range(1, 3):
    #     dataLists = ["results/ribo5/csv/sim100/modif{}/{}_merge.csv".format(modif, "_".join(i)) for i in drugNameList]
    #     saveName = "results/ribo5/images/sameDrug_sim100_modif{}.png".format(modif)
    #     createHeatmap(drugNameList, dataLists, [2, 2], saveName)


    ## differentDrug combination
    # drugNameList = list(itr.combinations(drugNames, 2))
    # for modif in range(1, 3):
    #     dataLists = ["results/ribo5/csv/sim100/modif{}/{}_merge.csv".format(modif, "_".join(i)) for i in drugNameList]
    #     saveName = "results/ribo5/images/diffDrug_sim100_modif{}.png".format(modif)
    #     createHeatmap(drugNameList, dataLists, [3, 2], saveName)

    # combinatorial simulation(virtual drug)
    ## merge DataFiles
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["30s", "30s"], ["30s", "50s"]]
    # nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for drugName in drugNameList for target in targetList]
    # for modif in range(1, 3):
    #     csvdir = "results/ribo5/csv/sim100_v/modif{}".format(modif)
    #     for name in nameList:
    #         dirName = "{}/{}".format(csvdir, "_".join(name))
    #         fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
    #         df = mergeResults(fileNameList)
    #         df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(name)), index=False)
            
    ## create Image
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["30s", "30s"], ["30s", "50s"]]
    # nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for target in targetList for drugName in drugNameList]
    # for modif in range(1, 3):
    #     dataLists = ["results/ribo5/csv/sim100_v/modif{}/{}_merge.csv".format(modif, "_".join(name)) for name in nameList]
    #     saveName = "results/ribo5/images/virtualDrug_sim100_modif{}.png".format(modif)
    #     createHeatmap(nameList, dataLists, [3, 2], saveName)

    # oldeval simulation
    ## merge DataFiles
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # csvdir = "results/ribo5/csv/old100"
    # for modif in range(1, 3):
    #     for drugName in drugNameList:
    #         dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(drugName))
    #         fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
    #         df = mergeResults(fileNameList)
    #         df.to_csv("{}_merge.csv".format(dirName), index=False)

    ## SameDrug combination
    # drugNameList = [[name, name] for name in drugNames]
    # csvdir = "results/ribo5/csv/old100"
    # for modif in range(1, 3):
    #     dataLists = ["{}/modif{}/{}_merge.csv".format(csvdir, modif, "_".join(i)) for i in drugNameList]
    #     saveName = "results/ribo5/images/sameDrug_old100_modif{}.png".format(modif)
    #     createHeatmap(drugNameList, dataLists, [2, 2], saveName, simType="epsilon")
    
    ## diffDrug combination
    # drugNameList = list(itr.combinations(drugNames, 2))
    # csvdir = "results/ribo5/csv/old100"
    # for modif in range(1, 3):
    #     dataLists = ["{}/modif{}/{}_merge.csv".format(csvdir, modif, "_".join(i)) for i in drugNameList]
    #     saveName = "results/ribo5/images/diffDrug_old100_modif{}.png".format(modif)
    #     createHeatmap(drugNameList, dataLists, [3, 2], saveName, simType="epsilon")

    # oldeval simulation (virtual drug)
    ## merge DataFiles
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["30s", "30s"], ["30s", "50s"]]
    # nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for drugName in drugNameList for target in targetList]
    # csvdir = "results/ribo5/csv/old100_v"
    # for modif in range(1, 3):
    #     for name in nameList:
    #         dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(name))
    #         fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
    #         df = mergeResults(fileNameList)
    #         df.to_csv("{}_merge.csv".format(dirName), index=False)
            
    ## create Image
    # drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    # targetList = [["30s", "30s"], ["30s", "50s"]]
    # nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for target in targetList for drugName in drugNameList]
    # csvdir = "results/ribo5/csv/old100_v"
    # for modif in range(1, 3):
    #     dataLists = ["{}/modif{}/{}_merge.csv".format(csvdir, modif, "_".join(name)) for name in nameList]
    #     saveName = "results/ribo5/images/virtualDrug_old100_modif{}.png".format(modif)
    #     createHeatmap(nameList, dataLists, [3, 2], saveName, simType="epsilon")
    
    # neweval simulation
    ## merge DataFiles
    # drugNameList = list(itr.combinations_with_replacement(drugNames, 2))
    # csvdir = "results/ribo5/csv/new100"
    # for modif in range(1, 3):
    #     for drugName in drugNameList:
    #         dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(drugName))
    #         fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(5)]
    #         df = mergeResults(fileNameList)
    #         df.to_csv("{}_merge.csv".format(dirName), index=False)

    ## SameDrug combination
    # drugNameList = [[name, name] for name in drugNames]
    # csvdir = "results/ribo5/csv/new100"
    # for modif in range(1, 3):
    #     dataLists = ["{}/modif{}/{}_merge.csv".format(csvdir, modif, "_".join(i)) for i in drugNameList]
    #     saveName = "results/ribo5/images/sameDrug_new100_modif{}.png".format(modif)
    #     createHeatmap(drugNameList, dataLists, [2, 2], saveName, xy=["Slope", "MidPoint"], simType="LinerType")
    
    ## diffDrug combination
    # drugNameList = list(itr.combinations(drugNames, 2))
    # csvdir = "results/ribo5/csv/new100"
    # for modif in range(1, 3):
    #     dataLists = ["{}/modif{}/{}_merge.csv".format(csvdir, modif, "_".join(i)) for i in drugNameList]
    #     saveName = "results/ribo5/images/diffDrug_new100_modif{}.png".format(modif)
    #     createHeatmap(drugNameList, dataLists, [3, 2], saveName, xy=["Slope", "MidPoint"], simType="LinerType")

    # neweval simulation (virtual drug)
    ## merge DataFiles
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["30s", "30s"], ["30s", "50s"]]
    nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for drugName in drugNameList for target in targetList]
    csvdir = "results/ribo5/csv/new100_v"
    for modif in range(1, 3):
        for name in nameList:
            dirName = "{}/modif{}/{}".format(csvdir, modif, "_".join(name))
            fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(5)]
            df = mergeResults(fileNameList)
            df.to_csv("{}_merge.csv".format(dirName), index=False)
            
    ## create Image
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["30s", "30s"], ["30s", "50s"]]
    nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for target in targetList for drugName in drugNameList]
    csvdir = "results/ribo5/csv/new100_v"
    for modif in range(1, 3):
        dataLists = ["{}/modif{}/{}_merge.csv".format(csvdir, modif, "_".join(name)) for name in nameList]
        saveName = "results/ribo5/images/virtualDrug_new100_modif{}.png".format(modif)
        createHeatmap(nameList, dataLists, [3, 2], saveName, xy=["Slope", "MidPoint"], simType="LinerType")
       


