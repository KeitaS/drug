#coding: utf-8
from modules import *

if __name__ == "__main__":
    drugNames = ["Streptmycin", "Kanamycin", "Tetracycline", "Chloramphenicol"]
    
    # single simulation 
    # singleDataPath = "results/ribo8/csv/single"
    # imagesDirPath = "results/ribo8/images"
    # plt.figure()
    # for index, drugName in enumerate(drugNames):
    #     plt.subplot(2, 2, index + 1)
    #     df = pd.read_csv("{}/{}.csv".format(singleDataPath, drugName))
    #     plt.plot(list(df["dose"]), list(df["growthRate"]))
    #     plt.title(drugName, fontsize=12)
    #     plt.xlabel("Dose", fontsize=12)
    #     plt.ylabel("Growth Rate", fontsize=12)
    # plt.tight_layout()
    # plt.savefig("{}/single.png".format(imagesDirPath), dpi=300)

    # combinatorial simulation
    ## merge DataFiles
    # drugNameList = itr.combinations_with_replacement(drugNames, 2)
    # csvdir = "results/ribo8/csv/double/normal"
    # for drugName in drugNameList:
    #     dirName = "{}/{}".format(csvdir, "_".join(drugName))
    #     fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
    #     df = mergeResults(fileNameList)
    #     df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(drugName)), index=False)

    ## SameDrug combination
    # drugNameList = [[name, name] for name in drugNames]
    # dataLists = ["results/ribo8/csv/double/normal/{}_merge.csv".format("_".join(i)) for i in drugNameList]
    # saveName = "results/ribo8/images/sameDrug.png"
    # createHeatmap(drugNameList, dataLists, [2, 2], saveName, ["dose2", "dose1"], "growthRate")

    ## differentDrug combination
    # drugNameList = list(itr.combinations(drugNames, 2))
    # dataLists = ["results/ribo8/csv/double/normal/{}_merge.csv".format("_".join(i)) for i in drugNameList]
    # saveName = "results/ribo8/images/diffDrug.png"
    # createHeatmap(drugNameList, dataLists, [3, 2], saveName, ["dose2", "dose1"], "growthRate")

    # combinatorial simulation(virtual drug)
    ## merge DataFiles
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["A", "A"], ["A", "B"], ["A", "C"]]
    csvdir = "results/ribo8/csv/double/virtual"
    for drugName in drugNameList:
        for target in targetList:
            dirName = "{}/{}".format(csvdir, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))]))
            fileNameList = ["{}/{}.csv".format(dirName, num) for num in range(101)]
            df = mergeResults(fileNameList)
            df.to_csv("{}/{}_merge.csv".format(csvdir, "_".join(["{}{}".format(drugName[i], target[i]) for i in range(len(drugName))])), index=False)

    ## create Image
    drugNameList = [["Streptmycin", "Streptmycin"], ["Streptmycin", "Chloramphenicol"], ["Chloramphenicol", "Chloramphenicol"]]
    targetList = [["A", "A"], ["A", "B"], ["A", "C"]]
    nameList = [["{}{}".format(drugName[0], target[0]), "{}{}".format(drugName[1], target[1])] for target in targetList for drugName in drugNameList]
    dataLists = ["results/ribo8/csv/double/virtual/{}_merge.csv".format("_".join(name)) for name in nameList]
    saveName = "results/ribo8/images/virtualDrug.png"
    createHeatmap(nameList, dataLists, [3, 3], saveName, xy=["dose2", "dose1"], simType="growth")

