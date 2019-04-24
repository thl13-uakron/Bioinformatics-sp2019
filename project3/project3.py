# data structures
## gene probes
def Probe(i="", d=""):
    return {"id":i, "desc":d, "expVals":{"AML":{}, "ALL":{}}}

# helper functions
## obtain list of genes for data file
def getProbeList(datasetFilename, classVectorFilename):
    # helper to clean up data lists
    def removeAll(l, c):
        for i in l:
            if i == c:
                l.remove(i)
          
    probes = []
    
    # open class vector file
    classVector = open(classVectorFilename)

    # class vector file format
    ## number of microarray tests (deliminator: ' ' after)
    ## number of (' ' after)
    ## number of (\n after)
    ## classification of each test sample (0 - ALL, 1 - AML) (' ' after)

    # get index data
    indices = classVector.readline()
    indices = indices.split(" ")
    testCount = int(indices[0])
    
    testClasses = classVector.readline()
    testClasses = testClasses.split(" ")
    testClasses.remove("\n")
    
    # open dataset file
    dataset = open(datasetFilename)
    
    # dataset file format:
    ## header fields (deliminator: \n after)
    ### description, accessor id (\t between)
    ### microarray test samples (\t\t between)
    ## list of test sample ids (\t\t between, \n after)
    ## number of genes (\n after)
    ## list of genes (\n between)
    ### data for header fields (\t between)
    ### test scores
    #### expression value (\t after)
    #### expression label (P, M, A) (\t after)
    #### next expression value, etc.
    
    # record header fields
    titleList = dataset.readline()
    titleList = titleList.split("\t")
    titleList.remove("\n")
    
    testList = dataset.readline()
    testList = testList.split("\t")
    testList.remove("\n")
    removeAll(testList, "")
    
    # record gene data
    probeCount = dataset.readline()
    probeStrings = dataset.readlines()
    # convert to objects, add to probe list
    for p in probeStrings:
        p = p.split("\t") # split up data elements
        p.remove("\n")

        # assign data to class members
        newProbe = Probe(p[1], p[0])
        
        pIndex = 2 # iterate through test values in data list
        while pIndex < len(p):
            testResult = {"score":p[pIndex], "label":p[pIndex + 1]}
            
            testIndex = int((pIndex - 2) / 2) # corresponding test sample
            testId = testList[testIndex] # get test name and classification
            if testClasses[testIndex] == '1': 
                newProbe["expVals"]["AML"][testId] = testResult
            else:
                newProbe["expVals"]["ALL"][testId] = testResult
            
            pIndex += 2
        
        probes.append(newProbe)

    # close files
    classVector.close()
    dataset.close()

    return probes

## write gene list to file in same format as original
def writeDatasetFile(probeList, newDatasetFilename, oldDatasetFilename, oldClassVectorFilename):
    return

# initial file data
trainingDatasetFile = "ALL_vs_AML_train_set_38_sorted.res"
trainingClassVector = "ALL_vs_AML_train_set_38_sorted.cls"
preprocessedTrainingDataset = "ALL_vs_AML_train_set_38_preprocessed.res"
geneList = getProbeList(trainingDatasetFile, trainingClassVector)

# preprocessing
threshold = 20
for gene in geneList:
    ## remove endrogenous control (housekeeping) genes
    if "endogenous control" in gene["desc"]:
        geneList.remove(gene)
    else:
        ## remove genes labelled "A" in all experiments
        ## replace low expression values with score threshold
        ## remove genes where max value is less than double min value
        toRemove = True
        minVal = threshold
        maxVal = threshold
        
        for testClass in gene["expVals"]:
            for testId in gene["expVals"][testClass]:
                score = gene["expVals"][testClass][testId]["score"]
                if int(score) < threshold:
                    score = threshold
                    gene["expVals"][testClass][testId]["score"] = score
                elif int(score) > maxVal:
                    maxVal = int(score)
                    
                if gene["expVals"][testClass][testId]["label"] != "A":
                    toRemove = False

        if maxVal < minVal * 2:
            toRemove = True
        if toRemove == True:
            geneList.remove(gene)
                    
## send preprocessed data to file
writeDatasetFile(geneList, preprocessedTrainingDataset, trainingDatasetFile, trainingClassVector)
