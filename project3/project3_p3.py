# data structures
from project3_p2 import Probe
from project3_p2 import getProbeList

# helper functions
## get manhattan distance between two genes
def getDistance(g1, g2):
    d = 0
    return d

# file data
trainingDatasetFile = ""
trainingClassVector = "ALL_vs_AML_train_set_38_sorted.cls"
trainingGeneList = getProbeList(trainingDatasetFile, trainingClassVector)

testingDatasetFile = "Leuk_ALL_AML.test.res"
testingClassVector = "Leuk_ALL_AML.test.cls"
testingGeneList = getProbeList(testingDatasetFile, testingClassVector)

# classification algorithm
## start with genes in training set
## take gene from testing set
## find genes in training set closest to selected gene
## give selected gene same classification as majority of nearby genes
## move gene to training set, repeat with rest of testing set
