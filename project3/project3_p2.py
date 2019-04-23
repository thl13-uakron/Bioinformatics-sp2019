# data structures
from project3 import Probe

# file data
trainingDatasetFile = ""
trainingClassVector = "ALL_vs_AML_train_set_38_sorted.cls"

# feature selection
## determine statistical difference in average expression value of each set
## sort genes by p-value from two-sample T test, select and save 50 lowest
