# data structures
## gene probes
class Probe:
    identifier = ""
    description = ""
    expressionVals = {};
    def __init__(i="", d="", expVals = {"AML":[], "ALL":[]}):
        identifier = i
        description = d
        expressionValues = expVals
    #

# helper functions
## obtain list of genes for data file
def getProbeList(datasetFilename, classVectorFilename):
    probes = []
    
    # open class vector file
    classVector = open(classVectorFilename)
    # get indicies
    
    # open dataset file
    dataset = open(datasetFilename)
    
    # record header fields as strings
    titleListString = dataset.readline()
    testListString = dataset.readline()
    # convert to lists
    titleList = titleListString.split("\t")
    testList = testListString.split("\t")
    
    # record gene data as string
    probeString = dataset.readlines()
    # convert to objects, add to list
    
    # close files
    classVector.close()
    dataset.close()

    return probes

## write gene list to file in same format as original

# initial file data
trainingDatasetFile = "ALL_vs_AML_train_set_38_sorted.res"
trainingClassVector = "ALL_vs_AML_train_set_38_sorted.cls"
geneList = getProbeList(trainingDataset, trainingClassVector)

# preprocessing
## remove endrogenous control (housekeeping) genes
## remove genes labelled "A" in all experiments
## replace low expression values with score threshold
## remove genes where max value is less than double min value
## send preprocessed data to file
