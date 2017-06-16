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
resolutionLevel='8k'
print('working with resolution level %s...'%resolutionLevel)
entropyThresholds={}
entropyThresholds['200']=1.76149386213
entropyThresholds['2k']=1.15216894394
entropyThresholds['4k']=0.837034672379  
entropyThresholds['8k']=0.452935086356  
entropyThresholds['16k']=0.0459240474862 
entropyThresholds['23k']=0.

dataFile='../../../data/case.1/original/GSE72056_melanoma_single_cell_revised_v2.txt'

immuneCellsDataFileLearning='../../../data/case.1/formatted/nonMalignant.%sgenes.data.learning.csv'%resolutionLevel
immuneCellsMetadataFileLearning='../../../data/case.1/formatted/nonMalignant.%sgenes.immuneMetadata.learning.csv'%resolutionLevel

immuneCellsDataFileTesting='../../../data/case.1/formatted/nonMalignant.%sgenes.data.prediction.csv'%resolutionLevel
tumorCellsMetadataFileTesting='../../../data/case.1/formatted/nonMalignant.%sgenes.tumorMetadata.prediction.csv'%resolutionLevel
immuneCellsMetadataFileTesting='../../../data/case.1/formatted/nonMalignant.%sgenes.immuneMetadata.prediction.csv'%resolutionLevel

testingNumCells=5000 # there is a total of 4645 cells, so if this number is larger than that, all cells will be included
testingNumVar=25000 # there is a total of 23686 genes, so if this number is larger than that, all genes will be included

learningTumors=['53','58','60','72','74','84','94']
testingTumors=['79','88','78','81','80','89']
selectedTumors=learningTumors+testingTumors

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

# 0.2. defining some variables
allGeneNames=[]
entropies=[]

# 0.3. removing previous data
if os.path.exists(immuneCellsDataFileLearning) == True:
    os.remove(immuneCellsDataFileLearning)
    os.remove(immuneCellsDataFileTesting)

# 0.4. starting metadata files
gMetaTumor=open(tumorCellsMetadataFileTesting,'w')
gMetaTumor.write('cell.instance,tumor.label\n')

gMetaImmune=open(immuneCellsMetadataFileLearning,'w')
gMetaImmune.write('cell.instance,immune.label\n')

gMetaImmuneTP=open(immuneCellsMetadataFileTesting,'w')
gMetaImmuneTP.write('cell.instance,immune.label\n')

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

    gIL=open(immuneCellsDataFileLearning,'a')
    gIT=open(immuneCellsDataFileTesting,'a')

    # adding cell labels
    for i in range(len(cellIDs)):
        if malignantLabels[i] == '1' and tumorLabels[i] in learningTumors and immuneLabels[i] != '0':
            gIL.write(',')
            gIL.write(cellIDs[i])
        if malignantLabels[i] == '1' and tumorLabels[i] in testingTumors and immuneLabels[i] != '0':
            gIT.write(',')
            gIT.write(cellIDs[i])
    gIL.write('\n')
    gIT.write('\n')

    # adding tumor cell labels
    for i in range(len(cellIDs)):
        if malignantLabels[i] == '1' and tumorLabels[i] in testingTumors and immuneLabels[i] != '0': 
            string2Write=cellIDs[i]+','+'Mel'+tumorLabels[i]+'\n'
            gMetaTumor.write(string2Write)
            selectedTumorsRanks[tumorLabels[i]]=selectedTumorsRanks[tumorLabels[i]]+1
    
    # adding immune cell type
    for i in range(len(cellIDs)):
        if malignantLabels[i] == '1' and tumorLabels[i] in learningTumors and immuneLabels[i] != '0':
            string2Write=cellIDs[i]+','+immuneCode[immuneLabels[i]]+'\n'
            gMetaImmune.write(string2Write)

        if malignantLabels[i] == '1' and tumorLabels[i] in testingTumors and immuneLabels[i] != '0':
            string2Write=cellIDs[i]+','+immuneCode[immuneLabels[i]]+'\n'
            gMetaImmuneTP.write(string2Write)

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
            gIL.write(geneName)
            gIT.write(geneName)
            for i in range(len(cellIDs)):
                if malignantLabels[i] == '1' and tumorLabels[i] in learningTumors and immuneLabels[i] != '0':
                    gIL.write(',')
                    gIL.write(expression[i])
                if malignantLabels[i] == '1' and tumorLabels[i] in testingTumors and immuneLabels[i] != '0':
                    gIT.write(',')
                    gIT.write(expression[i])

            gIL.write('\n')
            gIT.write('\n')
            lineCount=lineCount+1

gIL.close()
gIT.close()

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
