# data structures
from project3_p2 import Probe
from project3_p2 import getProbeList

# helper functions
## get list of samples and corresponding expression values from gene list
def getSampleList(geneList, restrictedIds=[]):
    sampleList = {}
    for gene in geneList:
        if gene["id"] not in restrictedIds:
            for testClass in gene[expVals]:
                # essentially "rotating" the score table portion
                # {geneId1:{sampleId1:score11, sampleId2, score12, etc.
                # -> {sampleId1:{geneId1:score11, geneId2:score21, etc.
                scoreList = gene[expVals][testClass]
                for sample in scoreList:
                    sampleList[sample][gene["id"]] = int(scoreList[sample]["score"])
    return sampleList

## get classifications of corresponding samples
def getSampleClasses(classVectorFilename):
    classVector = open(classVectorFilename)
    
    # extract line containing classifications, convert to list
    classVector.readline()
    sampleClasses = classVector.readline()
    sampleClasses = sampleClasses.replace("\n", "")
    sampleClasses = sampleClasses.split(" ")
        
    classVector.close()
    return sampleClasses

## get manhattan distance between two samples based on expression values
def getDistance(s1, s2):
    d = 0
    return d

# file data
trainingDatasetFile = "ALL_vs_AML_train_set_38_processed.res"
trainingClassVector = "ALL_vs_AML_train_set_38_sorted.cls"
trainingGeneList = getProbeList(trainingDatasetFile, trainingClassVector)

testingDatasetFile = "Leuk_ALL_AML_processed.test.res"
testingClassVector = "Leuk_ALL_AML.test.cls"
testingGeneList = getProbeList(testingDatasetFile, testingClassVector)

# expression data only compared for genes present in both datasets (optimization)
ignoredGenesList = []
for gene in testingGeneList:
    if gene not in trainingGeneList:
        ignoredGenesList.append(gene["id"])

# sample ids and expression values
trainingSampleList = getSampleList(trainingGeneList)
testingSampleList = getSampleList(testingGeneList, ignoredGenesList)

# classification algorithm
## start with samples in training set
## take samples from testing set
## find samples in training set closest to selected samples
## give selected sample same classification as majority of nearby sample
## move sample to training set, repeat with rest of testing set

# display results
## check accuracy based on actual classifications listed in file
