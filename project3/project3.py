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

def getProbeList(datasetFilename, classVectorFilename):
    probes = []
    # open class vector file
    # get indicies
    # open dataset file
    # record header fields
    # record gene data
    return probes

# initial file data
trainingDataset = "ALL_vs_AML_train_set_38_sorted.res"
trainingClassVector = ""
geneList = getProbelist(trainingDataset, trainingClassVector)

# preprocessing
## remove endrogenous control (housekeeping) genes
## remove genes labelled "A" in all experiments
## replace low expression values with score threshold
## remove genes where max value is less than double min value
## send preprocessed data to file
