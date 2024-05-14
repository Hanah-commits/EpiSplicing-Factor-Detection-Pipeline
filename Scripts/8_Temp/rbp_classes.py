import pandas as pd
import json

panther = pd.read_csv('./pantherGeneList.txt', delimiter='\t', header=None)
panther.columns = ['id', 'rbp', 'meta', 'gene_name', 'function', 'species']

def find_overlap(dictionary):
    overlaps = {}
    keys = list(dictionary.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            key1 = keys[i]
            key2 = keys[j]
            overlap = set(dictionary[key1]) & set(dictionary[key2])
            if overlap:
                overlaps[(key1, key2)] = list(overlap)
    return overlaps


## Functional classes: split into Splicing factors, RNA processing proteins, RNA METABOLISM proteins, misc, unknown functions
classes = {}
i = 0
for type in list(set(panther.function.values)):

    classes[type] = list(set(panther[panther.function == type].rbp.values))

classes['Unknown'] = list(set(panther[panther.function.isna()].rbp.values))


overlap = find_overlap(classes)

## Remove overlaps
for key_pair, common_values in overlap.items():
    key1, key2 = key_pair
    # print(f"Overlap between {key1} and {key2}: {common_values}")

    if key1 in ['RNA processing factor(PC00147)', 'RNA splicing factor(PC00148)', 'RNA metabolism protein(PC00031)']:
        classes[key2].remove(common_values[0])
        
    elif key2 in ['RNA processing factor(PC00147)', 'RNA splicing factor(PC00148)', 'RNA metabolism protein(PC00031)']:
        classes[key1].remove(common_values[0])

    elif key1 == 'Unknown':
        classes[key1].remove(common_values[0])
    elif key2 == 'Unknown':
        classes[key2].remove(common_values[0])

## Filter out keys with empty lists
classes = {key: value for key, value in classes.items() if value}

## Get Misc classes
misc = [k for k in classes.keys() if k not in ['RNA processing factor(PC00147)', 'RNA splicing factor(PC00148)', 'RNA metabolism protein(PC00031)', 'Unknown']]
classes['Misc'] = [value for key, values in classes.items() if key in misc for value in values] # Collapse values of specified keys

# Filter out keys in misc
classes = {key: value for key, value in classes.items() if key not in misc}

## write to JSON 
with open('HelperFunctions/RBP_Classes.json', 'w') as json_file:
    json.dump(classes, json_file, indent=4)  

