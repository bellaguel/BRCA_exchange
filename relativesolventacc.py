#!/usr/bin/env python3
#File: relativesolventacc.py
#Purpose: Create a pruned table with relevant information to detect clinical significance in relation to RSA in BRCA variants
#Author: Bella Guel
#Collaborators: Letitia Mueller 


import pandas as pd
import argparse
#import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def argParser():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--outputfile','-o',type=str,help='specify filename of your outputfile')
    return parser.parse_args()


args = argParser()
outfile = args.outputfile

rsa_adj ={
    'A':121,
    'R':265,
    'N':187,
    'D':187,
    'C':148,
    'E':214,
    'Q':214,
    'G':97,
    'H':216,
    'I':195,
    'L':191,
    'K':230,
    'M':203,
    'F':228,
    'P':154,
    'S':143,
    'T':163,
    'W':264,
    'Y':255,
    'V':165
} ## Dictionary used to calculate relative solvent accessibility from ACC value from DSSP

BRCA_info = pd.read_csv('~/Desktop/SIP/BRCA/brca_info.csv', header=0) #reads in CSV file into pandas dataframe, calls datafram BRCA_info

BRCA_info = BRCA_info[BRCA_info.mupit_structure != '-'] # taking all mupit structures that are not equal to a hypen
BRCA_info = BRCA_info[BRCA_info.mupit_structure != 'fENSP00000380152_7'] #do not want this varrient
BRCA_info = BRCA_info[BRCA_info.mupit_structure != ''] #dropping no info
BRCA_info = BRCA_info.dropna(subset=['Pos']) # dropping Not applicable values from pos column
BRCA_info = BRCA_info[BRCA_info.Pos != '-'] #taking all positive values that are not equal to a hypen


combined_info = pd.DataFrame(columns=['amino acid', 'protein structure pos', 'genomic coord','mupit pos','relativesolventacc','protein_prior','clinicalsig', 'PPclinsig', 'grouping'])
# create a new pandas dataframe to name columns of final table

combined_info['clinicalsig']=BRCA_info['Clinical_Significance_ClinVar'] # move column "clinical sig clinvar" from BRCA_info dataframe to combined_info dataframe under column name "clinicalsig"
combined_info['protein_prior']=BRCA_info['proteinPrior'] #moving column "proteinPrior" from BRCA info to combined_info "protein prior"
combined_info['genomic coord']=BRCA_info['Pos'] # moving "Pos" column from BRCA_info to "genomic coord" in combined_info


mupit_4igk = pd.read_csv('~/Desktop/SIP/BRCA/4igk.csv',header=0) # reading in the mupit csv files, and ignoring header
mupit_1t15 = pd.read_csv('~/Desktop/SIP/BRCA/1t15.csv',header=0)
mupit_1jm7 = pd.read_csv('~/Desktop/SIP/BRCA/1jm7.csv',header=0)


dssp_4igk = pd.read_csv('~/Desktop/SIP/BRCA/4igk_edit.dssp', sep='\s+', header=0) # reading in the dssp csv files, denote space separator (not comma separator), ignoring header names from priorr chart by denoting 0
dssp_1t15 = pd.read_csv('~/Desktop/SIP/BRCA/1t15_edit.dssp',sep='\s+',header=0)
dssp_1jm7 = pd.read_csv('~/Desktop/SIP/BRCA/1jm7_edit.dssp',sep='\s+',header=0)


helper_4igk = pd.DataFrame(columns=['amino acid','relativesolventacc','protein structure pos','pos1','pos2','pos3']) # creating pandas dataframe as a helper to combine dssp data and mupit data
helper_1t15 = pd.DataFrame(columns=['amino acid','relativesolventacc','protein structure pos','pos1','pos2','pos3'])
helper_1jm7 = pd.DataFrame(columns=['amino acid','relativesolventacc','protein structure pos','pos1','pos2','pos3'])

helpers = [helper_4igk,helper_1t15,helper_1jm7]
dsspfiles = [dssp_4igk,dssp_1t15,dssp_1jm7]
mupitfiles = [mupit_4igk,mupit_1t15,mupit_1jm7]

#
# for help in helpers:  ###### not sure why i cant put this into a for loop, if you figure it out youre a freaking superstar
#     for dsspfile in dsspfiles:
#         help['amino acid']=dsspfile['AA']
#         help['relativesolventacc'
#         ]
## I have tried to shorten this but it results in a key error


## combining mupit data and dssp data into helper dataframe for each protein structure
print("adding columns from mupit and dssp 1") # relocating info from helper files into main data frame equivalents
helper_4igk['amino acid'] = dssp_4igk['AA']
helper_4igk['relativesolventacc'] = dssp_4igk['ACC']
helper_4igk['protein structure pos'] = mupit_4igk['seqRes']
helper_4igk['pos1'] = mupit_4igk['pos1']
helper_4igk['pos2'] = mupit_4igk['pos2']
helper_4igk['pos3'] = mupit_4igk['pos3']

print("adding columns from mupit and dssp 2")

helper_1t15['amino acid'] = dssp_1t15['AA'] #relocating info from helper files into main data frame equivalents
helper_1t15['relativesolventacc'] = dssp_1t15['ACC']
helper_1t15['protein structure pos'] = mupit_1t15['seqRes']
helper_1t15['pos1'] = mupit_1t15['pos1']
helper_1t15['pos2'] = mupit_1t15['pos2']
helper_1t15['pos3'] = mupit_1t15['pos3']

print("adding columns from mupit and dssp 3") #relocating info from helper files into main data frame equivalents
helper_1jm7['amino acid'] = dssp_1jm7['AA']
helper_1jm7['relativesolventacc'] = dssp_1jm7['ACC']
helper_1jm7['protein structure pos'] = mupit_1jm7['seqRes']
helper_1jm7['pos1'] = mupit_1jm7['pos1']
helper_1jm7['pos2'] = mupit_1jm7['pos2']
helper_1jm7['pos3'] = mupit_1jm7['pos3']



print("matching genomic coordinates ") # stratigic print statement to see when data starts to match in main dataframe

# matching genomic coordinate positions from helper dataframe to combined_info dataframe --> working in rows
for filename in [helper_4igk,helper_1jm7,helper_1t15]:
    for index,row in filename.iterrows():
        for i,r in combined_info.iterrows():
            if int(r['genomic coord'])==row['pos1']:
                r['amino acid']=row['amino acid']
                r['protein structure pos']= row['protein structure pos']
                r['relativesolventacc']=row['relativesolventacc']
            if int(r['genomic coord'])==row['pos2']:
                r['amino acid']=row['amino acid']
                r['protein structure pos']= row['protein structure pos']
                r['relativesolventacc']=row['relativesolventacc']
            if int(r['genomic coord'])==row['pos3']:
                r['amino acid']=row['amino acid']
                r['protein structure pos']= row['protein structure pos']
                r['relativesolventacc']=row['relativesolventacc']
    print("matching genomic coordinates 1")


print("cleaning up final data")
combined_info = combined_info.dropna(subset=['relativesolventacc']) #removing missing values within relativesolventacc

combined_info= combined_info[combined_info.clinicalsig != '-'] #repeating process of mupit structres in order to get desired data
combined_info= combined_info[combined_info.clinicalsig != 'Uncertain_significance']
combined_info= combined_info[combined_info.clinicalsig != 'not_provided']
combined_info= combined_info[combined_info.relativesolventacc != '']

# combined_info[['clinicalsig','sig1','sig2','sig3']] = combined_info['clinicalsig'].str.split(',',expand=True)

print("calculating relative solvent accessibility")
for i,r in combined_info.iterrows():
    r['relativesolventacc']= r['relativesolventacc']/rsa_adj[r['amino acid']] # calculating relative solvent accessibility based off baseline


print("separating clinvar clinical significance")
combined_info1 = combined_info['clinicalsig'].apply(lambda x: pd.Series(x.split(','))) # splitting on a comma
#lambda --> function that does not require a return to make it split on the comma
combined_info1 = combined_info1.rename(index=str, columns={0: "sig0", 1: "sig1", 2: "sig2", 3:"sig3", 4: "sig4"}) #adding new columns that come with the split


# combined_info = pd.concat([combined_info, combined_info1], axis=1) #linking up sigs in the combinded info
print("creating separate columns")
#
# for col_name in [sig0,sig1,sig2,sig3,sig4]:
#     for string_col_name in ['sig0','sig1','sig2','sig3','sig4']:

combined_info = combined_info.assign(sig0=combined_info1['sig0'].values) #combined_info from splits with other data
    # combined_info[string_col_name].replace('', np.nan, inplace=True)
combined_info = combined_info.assign(sig1=combined_info1['sig1'].values)
combined_info = combined_info.assign(sig2=combined_info1['sig2'].values)
combined_info = combined_info.assign(sig3=combined_info1['sig3'].values)
combined_info = combined_info.assign(sig4=combined_info1['sig4'].values)


# with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
#     print(combined_info)
print("replacing NA values")

# for string_col_name in
combined_info['sig0'].replace('', np.nan, inplace=True)
combined_info['sig1'].replace('', np.nan, inplace=True)
combined_info['sig2'].replace('', np.nan, inplace=True)
combined_info['sig3'].replace('', np.nan, inplace=True)
combined_info['sig4'].replace('', np.nan, inplace=True)

#calling N/A values that we do not need and replacing with empty strings
#making searching through data to classify easier

print("classifying clinical significance")
for index,row in combined_info.iterrows():
    cat_list = [] #making empty list in order
    cat_list.append(row['sig0'])
    cat_list.append(row['sig1'])
    cat_list.append(row['sig2'])
    cat_list.append(row['sig3'])
    cat_list.append(row['sig4'])
    # for i in range(0,5):
    #     print(row[i])
    #     cat_list.append(row[i])
    # print(cat_list)
    val = '0'
    if 'Benign' in cat_list:
        val = 'Benign'
        # print("assigned benign")
    elif 'Likely_benign' in cat_list: #assigning values for dataframe based on combinded info
        val = 'Benign'
    elif 'Pathogenic' in cat_list:
        val = 'Pathogenic'
        # print("assigned pathogenic")
    elif 'Likely_pathogenic' in cat_list:
        val = 'Pathogenic'
        # print("assigned pathogenic")

    else:
        val = "-"
        # print("assigned -")
    # print(val)
    row['clinicalsig']=val
combined_info = combined_info[combined_info.clinicalsig != "-"]

print("dumping to csv")


print("making pp clingsig")

combined_info = combined_info[combined_info.protein_prior != "-"]

for index, row in combined_info.iterrows():
    if float(row['protein_prior']) <0.30:
        print("if 1 working")
        row['PPclinsig'] = 'Benign'
    else:
        print("if 2 working")
        row['PPclinsig'] = 'Pathogenic'

for index, row in combined_info.iterrows():
    if row['PPclinsig'] and row['clinicalsig'] == 'Pathogenic':
        row['grouping'] = 'true positive'

for index, row in combined_info.iterrows():
    if row['PPclinsig'] and row['clinicalsig'] == 'Benign':
        row['grouping'] = 'true benign'

for index, row in combined_info.iterrows():
    if row['PPclinsig']== "Benign" and row['clinicalsig'] == 'Pathogenic':
        row['grouping'] = 'false negative'

for index, row in combined_info.iterrows():
    if row['PPclinsig'] == "Pathogenic" and row['clinicalsig'] == 'Benign':
        row['grouping'] = 'false positive'

combined_info.to_csv('table.csv')
print(combined_info)

BRCA_info = BRCA_info.dropna(subset=['Pos']) #removing values within subest

# fig1, ax1 = plt.subplots()
# ax1.set_title('Benign and Pathogenic Relative Solvent Accessibility (RSA)')
# ax1.boxplot(combined_info["relativesolventacc"])
# fig1

bp= combined_info.boxplot(by= "clinicalsig", column= ["relativesolventacc"], grid= True)
bp.get_figure().suptitle('Benign Varients Have Higher Relative Solvent Accessiblity (RSA)')
plt.savefig('Guel_Bella_boxplotnew.png',dpi=600)

bp= combined_info.boxplot(by= "grouping", column= ["relativesolventacc"], grid= True)
bp.get_figure().suptitle('RSA grouped by clinical significance predicted by protien priorr')
plt.savefig('Guel_Bella_protienpriorr_vs_RSA.png',dpi=600)
