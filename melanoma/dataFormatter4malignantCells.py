################################################################
# @Course: Systems Biology of Disease                           #
# @Python Script: dataFormatter4malignantCells.py               #
# @Author: Adrian de Lomana                                     #
# @Sponsored by:                                                #
# Institute for Systems Biology                                 #
# 401 Terry Ave N                                               #
# Seattle, WA 98109                                             #
#                                                               #
# Copyright (C) 2016 by Institute for Systems Biology,          #
# Seattle, Washington, USA.  All rights reserved.               #
#                                                               #
# This source code is distributed under GNU GPL v3.0            #
# the text of which is available at:                            #
# http://choosealicense.com/licenses/gpl-3.0/                   #
#################################################################

###
### this script performs takes the original data from PMID 27124452 and formats the data appropriately for analysis
###

import sys,os,scipy,math,numpy
import scipy.stats

def entropyCalculator(v):

    # calculating the probability distribution
    k=numpy.arange(0.,16.+1.)
    n,bins=numpy.histogram(v,bins=k)

    y=[]
    y=numpy.array(n)
    y=y/float(sum(y))

    s=scipy.stats.entropy(y)

    return s

# 0. user defined variables and paths
dataFile='data/original/GSE72056_melanoma_single_cell_revised_v2.txt'
tumorCellsFile='data/formatted/tumorCells.data.csv'
tumorCellsMetadataFile='data/formatted/tumorCells.metadata.csv'

selectedTumors=['79','88','84','78','81','80']

testingNumCells=5000 # there are 4645 cells in total
testingNumVar=24000   # there are 23k genes
# 5k and 16k is good

entropyThreshold=0.0692093239663
#entropyThreshold=0.

repeatedGenes=['MARCH1','MARCH2']

if os.path.exists(tumorCellsFile) == True:
    os.remove(tumorCellsFile)

gMetaTumor=open(tumorCellsMetadataFile,'w')
gMetaTumor.write('cell.instance,tumor.label\n')

allGeneNames=[]

# 1. slitting expression
with open(dataFile, 'r') as f:

    # dealing with headers
    cellIDs=f.readline().split('\t')
    cellIDs[-1]=cellIDs[-1].replace('\n','')
    del cellIDs[0]

    tumorLabels=f.readline().split('\t')
    tumorLabels[-1]=tumorLabels[-1].replace('\n','')
    del tumorLabels[0]

    malignantLabels=f.readline().split('\t')
    malignantLabels[-1]=malignantLabels[-1].replace('\n','')
    del malignantLabels[0]

    ### TBR
    malignantLabels=malignantLabels[:testingNumCells]
    tumorLabels=tumorLabels[:testingNumCells]
    ####
    
    gT=open(tumorCellsFile,'a')

    # adding cell labels
    for i in range(len(malignantLabels)):
        if malignantLabels[i] == '2' and tumorLabels[i] in selectedTumors:
            gT.write(',')
            gT.write(cellIDs[i])
    gT.write('\n')

    # adding tumor cell labels
    for i in range(len(malignantLabels)):
        if malignantLabels[i] == '2' and tumorLabels[i] in selectedTumors: 
            string2Write=cellIDs[i]+','+'Mel'+tumorLabels[i]+'\n'
            gMetaTumor.write(string2Write)
    
    ## # adding immune cell type
    ## gI.write('immuneLabel,')
    ## for i in range(len(malignantLabels)):
    ##     if malignantLabels[i] == '1':
    ##         gI.write(immuneLabels[i])
    ##         gI.write(',')
    ## gI.write('\n')

    # 2. dealing with data


    ### TBR
    lineCount=0
    ###

    entropies=[]
    for line in f:
        vector=line.split('\t')
        geneName=vector[0]
        expression=vector[1:]
        expression[-1]=expression[-1].replace('\n','')

        expressionValues=[float(element) for element in expression]
        s=entropyCalculator(expressionValues)
        entropies.append(s)

        if geneName not in allGeneNames:
            allGeneNames.append(geneName)
        else:
            print(geneName)
            print('warning, geneName found twice...')

        ### TBR
        #expression=expression[:testingNumCells]
        #if lineCount > testingNumVar:
        #    break
        ###
        
        # writing the gene name
        if s > entropyThreshold and geneName not in repeatedGenes:
            gT.write(geneName)

            for i in range(len(expression)):
                if malignantLabels[i] == '2' and tumorLabels[i] in selectedTumors: 
                    gT.write(',')
                    gT.write(expression[i])

            gT.write('\n')

gT.close()

entropies.sort(reverse=True)
threshold=entropies[16000] # 0.0692093239663
print(threshold)
