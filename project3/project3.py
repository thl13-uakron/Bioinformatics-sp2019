# data structures
## gene probes
class Probe:
    expressionVals = {"AML":[], "ALL":[]};
    #

# initial file data
filename = "ALL_vs_AML_train_set_38_sorted.res"
geneList = []

# preprocessing
## remove endrogenous control (housekeeping) genes
## remove genes labelled "A" in all experiments
## replace low expression values with score threshold
## remove genes where max value is less than double min value
## send preprocessed data to file
