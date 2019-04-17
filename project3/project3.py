# data structures

## gene probes
class Probe:
    expressionVals = {"AML":[], "ALL":[]};
    #
 

# functions

## data preprocessing
### remove endrogenous control (housekeeping) genes
### remove genes labelled "A" in all experiments
### replace low expression values with score threshold
### remove genes where max value is less than double min value
### send preprocessed data to file

## feature selection
### determine statistical difference in average expression value of each set
### remove genes where p>0.05 after two-sample T test


# driver program

## process file data
