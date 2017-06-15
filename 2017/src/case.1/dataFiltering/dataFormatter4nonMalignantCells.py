#################################################################
# @Course: Systems Biology of Disease  2017                     #
# @Python Script: dataFormatter4malignantCells.py               #
# @Author: Adrian de Lomana                                     #
#                                                               #
# This source code is distributed under GNU GPL v3.0            #
# the text of which is available at:                            #
# http://choosealicense.com/licenses/gpl-3.0/                   #
#################################################################

###
### this script performs takes the original data from PMID 27124452 and formats the data appropriately for analysis for immune cells
###

import sys,os,scipy,math,numpy,codecs
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

# 0.1. user defined variables and paths
resolutionLevel='23k'
print('working with resolution level %s...'%resolutionLevel)
entropyThresholds={}
entropyThresholds['200']=1.76196273819
entropyThresholds['2k']=1.15709452349
entropyThresholds['4k']=0.841156501026  
entropyThresholds['8k']=0.456759266159  
entropyThresholds['16k']=0.0459240474862 
entropyThresholds['23k']=0.

dataFile='../../../data/case.1/original/GSE72056_melanoma_single_cell_revised_v2.txt'

immuneCellsFile='../../../data/case.1/formatted/nonMalignant.%sgenes.data.csv'%resolutionLevel
immuneCellsMetadataFile='../../../data/case.1/formatted/nonMalignant.%sgenes.immuneMetadata.csv'%resolutionLevel
tumorCellsMetadataFile='../../../data/case.1/formatted/nonMalignant.%sgenes.tumorMetadata.csv'%resolutionLevel
centroidsFileName='../../../data/case.1/formatted/nonMalignant.immuneTypes.median.%sgenes.csv'%resolutionLevel


testingNumCells=5000 # there is a total of 4645 cells, so if this number is larger than that, all cells will be included
testingNumVar=25000 # there is a total of 23686 genes, so if this number is larger than that, all genes will be included

selectedTumors=['53','58','60','72','74','79','80','84','88','89','94']
selectedTumorsRanks={}
for patient in selectedTumors:
    selectedTumorsRanks[patient]=0

immuneCode={}
immuneCode['0']='no class'
immuneCode['1']='T-cells'
immuneCode['2']='B-cells'
immuneCode['3']='Macrophages'
immuneCode['4']='Endo.'
immuneCode['5']='CAFs'
immuneCode['6']='NK'

signature={}
signature['1']={}
signature['2']={}
signature['3']={}
signature['4']={}
signature['5']={}
signature['6']={}

# 0.2. defining some variables
allGeneNames=[]
entropies=[]

# 0.3. removing previous data
if os.path.exists(immuneCellsFile) == True:
    os.remove(immuneCellsFile)

# 0.4. starting metadata files
gMetaTumor=open(tumorCellsMetadataFile,'w')
gMetaTumor.write('cell.instance,tumor.label\n')

gMetaImmune=open(immuneCellsMetadataFile,'w')
gMetaImmune.write('cell.instance,immune.label\n')

# 1. slitting expression
with open(dataFile, 'r') as f:

    # dealing with headers
    originalCellIDs=f.readline().split('\t')
    originalCellIDs[-1]=originalCellIDs[-1].replace('\n','')
    del originalCellIDs[0]
    cellIDs=['cell'+str(i+1) for i in range(len(originalCellIDs))]

    tumorLabels=f.readline().split('\t')
    tumorLabels[-1]=tumorLabels[-1].replace('\n','')
    del tumorLabels[0]

    malignantLabels=f.readline().split('\t')
    malignantLabels[-1]=malignantLabels[-1].replace('\n','')
    del malignantLabels[0]

    immuneLabels=f.readline().split('\t')
    immuneLabels[-1]=immuneLabels[-1].replace('\n','')
    del immuneLabels[0]

    ### dealing with trimmed versions of the data
    cellIDs=cellIDs[:testingNumCells]
    tumorLabels=tumorLabels[:testingNumCells]
    malignantLabels=malignantLabels[:testingNumCells]
    immuneLabels=immuneLabels[:testingNumCells]
    ####

    gI=open(immuneCellsFile,'a')

    # adding cell labels
    for i in range(len(cellIDs)):
        if malignantLabels[i] == '1' and tumorLabels[i] in selectedTumors and immuneLabels[i] != '0':
            gI.write(',')
            gI.write(cellIDs[i])
    gI.write('\n')

    # adding tumor cell labels
    for i in range(len(cellIDs)):
        if malignantLabels[i] == '1' and tumorLabels[i] in selectedTumors and immuneLabels[i] != '0': 
            string2Write=cellIDs[i]+','+'Mel'+tumorLabels[i]+'\n'
            gMetaTumor.write(string2Write)
            selectedTumorsRanks[tumorLabels[i]]=selectedTumorsRanks[tumorLabels[i]]+1
    
    # adding immune cell type
    for i in range(len(cellIDs)):
        if malignantLabels[i] == '1' and tumorLabels[i] in selectedTumors and immuneLabels[i] != '0':
            string2Write=cellIDs[i]+','+immuneCode[immuneLabels[i]]+'\n'
            gMetaImmune.write(string2Write)

    # 2. dealing with data

    ### dealing with trimmed versions of the data
    lineCount=0
    ###

    for line in f:
        vector=line.split('\t')
        geneName=vector[0]
        expression=vector[1:]
        expression[-1]=expression[-1].replace('\n','')

        trimmedExpression=[]
        for i in range(len(cellIDs)):
            if malignantLabels[i] == '1' and tumorLabels[i] in selectedTumors and immuneLabels[i] != '0':
                trimmedExpression.append(expression[i])
        expressionValues=[float(element) for element in trimmedExpression]
        s=entropyCalculator(expressionValues)
        entropies.append(s)
        duplicate=False

        if geneName not in allGeneNames:
            allGeneNames.append(geneName)
        else:
            print(geneName)
            print('warning, geneName found twice...')
            duplicate=True

        ### dealing with trimmed versions of the data
        if lineCount == testingNumVar:
            break
        ###
        
        # writing the gene name
        if duplicate == False and s > entropyThresholds[resolutionLevel]:
            gI.write(geneName)
            for i in range(len(cellIDs)):
                if malignantLabels[i] == '1' and tumorLabels[i] in selectedTumors and immuneLabels[i] != '0':
                    gI.write(',')
                    gI.write(expression[i])

                    if immuneLabels[i] != '0':
                        if geneName not in signature[immuneLabels[i]].keys():
                            signature[immuneLabels[i]][geneName]=[expression[i]]
                        else:
                            signature[immuneLabels[i]][geneName].append(expression[i])
                    
            gI.write('\n')
            lineCount=lineCount+1

gI.close()

# 2. some prints about entropy thresholds
print('printing entropies...')
entropies.sort(reverse=True)

size=200
threshold=entropies[size] 
print(size,threshold)

size=2000
threshold=entropies[size] 
print(size,threshold)

size=4000
threshold=entropies[size] 
print(size,threshold)

size=8000
threshold=entropies[size] 
print(size,threshold)

size=16000
threshold=entropies[size] 
print(size,threshold)

# 3. writing the median of each cell type
cellTypeCodes=signature.keys()
selectedGenes=signature['1'].keys()
g=open(centroidsFileName,'w')
for cellTypeCode in cellTypeCodes:
    g.write(',%s'%immuneCode[cellTypeCode])
g.write('\n')
for geneID in selectedGenes:
    g.write('%s'%geneID)
    for cellTypeCode in cellTypeCodes:
        values=signature[cellTypeCode][geneID]
        expression=[float(element) for element in values]
        medianValue=numpy.median(expression)
        string2Write=','+str(medianValue)
        g.write(string2Write)
    g.write('\n')
g.close()

print(selectedTumorsRanks)

