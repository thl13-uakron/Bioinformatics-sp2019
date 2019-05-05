import pandas as pd
pd.set_option("display.max_colwidth", 10000)

train_set = r"data_set_ALL_AML_train.tsv"
test_set = r"data_set_ALL_AML_independent.tsv"
train_vec = r"ALL_vs_AML_train_set_38_sorted.cls"
test_vec = r"Leuk_ALL_AML.test.cls"

train_df = pd.read_csv(train_set, '\t')
test_df = pd.read_csv(test_set, '\t')

# Eliminate the “endogenous control” genes (housekeeping genes)
num_control = 0
for index, row in train_df.iterrows():
    if 'endogenous control' in row.name:
        train_df.drop(index, inplace=True)
        num_control += 1
print('Num control dropped:', num_control)
    
train_df = train_df.fillna(20)

# Eliminate the genes with all As across the experiments
# Train 1 -> 38
expr_range = range(1, 39)
num_all_as = 0
for index, row in train_df.iterrows():
    all_as = True
    for i in expr_range:
        if row[str(i)] != 'A':
            all_as = False
            break
    if all_as == True:
        num_all_as += 1
        try:
            train_df.drop(index, inplace=True)
        except:
            pass
print('Num all A\'s dropped:', num_all_as)


# Replace all the expression values below some threshold cut-off value to that threshold value (pick 20 to be the threshold cut-off value)
thresh = 20
expr_range = range(0, 38)
num_threshed = 0
for index, row in train_df.iterrows():
    for i in expr_range:
        if i == 0:
            key='call'
        else:
            key = f'call.{i}'
        if row[key] < thresh:
            num_threshed += 1
            train_df.loc[index, key] = thresh
            
        
print('Num under threshold:', num_threshed)

#stat of aml
switch = 28

header = "ID"
for i in range(1, 39):
    header += f"\tExp{i} "
    if i >= switch:
        header += "(AML)"
    else:
        header += "(ALL)"

body = ''

for index, row in train_df.iterrows():
    entry = ""
    first = row['Gene Description']
    vals = []
    for i in expr_range:
        if i == 0:
            key='call'
        else:
            key = f'call.{i}'
        vals.append(row[key])
    entry += first
    for v in vals:
        entry += f"\t{v}"
    body = body + entry + '\n'
        

out = header + '\n' + body

with open('train_processed.txt', 'w') as outfile:
    outfile.write(out)