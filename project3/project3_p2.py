# data structures
from project3 import Probe
from project3 import getProbeList
from project3 import writeDatasetFile

# helper functions
## two-sample t test
def TTest(l1, l2):
    p = 0
    return p

# file data (training set)
trainingDatasetFile = "ALL_vs_AML_train_set_38_preprocessed.res"
trainingClassVector = "ALL_vs_AML_train_set_38_sorted.cls"
processedTrainingDataset = "ALL_vs_AML_train_set_38_processed.res"
trainingGeneList = getProbeList(trainingDatasetFile, trainingClassVector)

# feature selection
## determine statistical difference in average expression value of training set
for gene in trainingGeneList:
    AMLset = []
    ALLset = []
    for testId in gene["expVals"]["AML"]:
        AMLset.append(int(gene["expVals"]["AML"][testId]["score"]))
    for testId in gene["expVals"]["ALL"]:
        ALLset.append(int(gene["expVals"]["ALL"][testId]["score"]))
    gene["p-value"] = TTest(AMLset, ALLset)
    
## sort genes by p-value from two-sample T test, select and save 50 lowest to new file

# file data (testing set)
testingDatasetFile = "Leuk_ALL_AML.test.res"
testingClassVector = "Leuk_ALL_AML.test.cls"
processedTestingDataset = "Leuk_ALL_AML_processed.test.res"
testingGeneList = getProbeList(testingDatasetFile, testingClassVector)

# processing
## replace low values with score threshold
threshold = 20
for gene in testingGeneList:
    for testClass in gene["expVals"]:
        for testId in gene["expVals"][testClass]:
            if int(gene["expVals"][testClass][testId]["score"]) < threshold:
                gene["expVals"][testClass][testId]["score"] = str(threshold)
## save to new file
writeDatasetFile(testingGeneList, processedTestingDataset, testingDatasetFile)    
