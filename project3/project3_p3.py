# data structures
from project3_p2 import Probe
from project3_p2 import getProbeList

# helper functions
## get list of samples and corresponding expression values from gene list
def getSampleList(geneList, restrictedIds=[]):
    sampleList = {}
    for gene in geneList:
        if gene["id"] not in restrictedIds:
            for testClass in gene["expVals"]:
                # essentially "rotating" the score table portion
                # {geneId1:{sampleId1:score11, sampleId2, score12, etc.
                # -> {sampleId1:{geneId1:score11, geneId2:score21, etc.
                scoreList = gene["expVals"][testClass]
                for sample in scoreList:
                    if sample not in sampleList:
                        sampleList[sample] = {"score":{}}
                    sampleList[sample]["score"][gene["id"]] = int(scoreList[sample]["score"])
    return sampleList

## get classifications of corresponding samples
def getSampleClasses(datasetFilename, classVectorFilename):
    classVector = open(classVectorFilename)
    dataset = open(datasetFilename)
    
    # get sample ids
    dataset.readline()
    idList = dataset.readline()
    idList = idList.split("\t")
    # get classifications
    classVector.readline()
    classList = classVector.readline()
    classList = classList.split(" ")
    # pair corresponding data values
    sampleClasses = {}
    for i, c in zip(idList, classList):
        if i == "" or i == "\n":
            idList.remove(i)
        if c == "" or c == "\n":
            classList.remove(c)
    for i, c in zip(idList, classList):
        sampleClasses[i] = c
        
    classVector.close()
    dataset.close()
    return sampleClasses

## get manhattan distance between two samples based on expression values
def getDistance(s1, s2):
    d = 0
    for gene in s1["score"]:
        d += abs(s1["score"][gene] - s2["score"][gene])
    return d

# file data
trainingDatasetFile = "ALL_vs_AML_train_set_38_processed.res"
trainingClassVector = "ALL_vs_AML_train_set_38_sorted.cls"
trainingGeneList = getProbeList(trainingDatasetFile, trainingClassVector)

testingDatasetFile = "Leuk_ALL_AML_processed.test.res"
testingClassVector = "Leuk_ALL_AML.test.cls"
testingGeneList = getProbeList(testingDatasetFile, testingClassVector)

# expression data only compared for genes present in both datasets (optimization)
trainingGeneIds = []
for gene in trainingGeneList:
    trainingGeneIds.append(gene["id"])
ignoredGenesList = []
for gene in testingGeneList:
    if gene["id"] not in trainingGeneIds:
        ignoredGenesList.append(gene["id"])

# sample ids and expression values
trainingSampleList = getSampleList(trainingGeneList)
testingSampleList = getSampleList(testingGeneList, ignoredGenesList)

# sample classifications for training set
trainingSampleClasses = getSampleClasses(trainingDatasetFile, trainingClassVector)
for sample in trainingSampleClasses:
    trainingSampleList[sample]["class"] = trainingSampleClasses[sample]

# classification algorithm
## start with samples in training set
## take samples from testing set
## find samples in training set closest to selected samples
## give selected sample same classification as majority of nearby sample
## move sample to training set, repeat with rest of testing set
k = 3
for sample1 in testingSampleList:
    distances = {}
    nearestNeighbors = []
    for sample2 in trainingSampleList:
        distances[sample2] = getDistance(testingSampleList[sample1], trainingSampleList[sample2])
        nearestNeighbors.append(sample2)
    nearestNeighbors.sort(key=lambda s:distances[s])
    nearestNeighbors = nearestNeighbors[:3]

# display results
## check accuracy based on actual classifications listed in file
