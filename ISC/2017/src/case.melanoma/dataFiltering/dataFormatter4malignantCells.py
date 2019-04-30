#################################################################
# @Course: Systems Biology of Disease 2017                      #
# @Python Script: dataFormatter4malignantCells.py               #
# @Author: Adrian de Lomana                                     #
#                                                               #
# This source code is distributed under GNU GPL v3.0            #
# the text of which is available at:                            #
# http://choosealicense.com/licenses/gpl-3.0/                   #
#################################################################

###
### this script performs takes the original data from PMID 27124452 and formats the data appropriately for analysis of malignant cells
### original data file found at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056
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

# 0.1. user defined variables and paths
resolutionLevel='2k'
print('working with resolution level %s...'%resolutionLevel)
entropyThresholds={}
entropyThresholds['200']=1.76282876498
entropyThresholds['2k']=1.4685431404109979
entropyThresholds['4k']=1.17578271113  
entropyThresholds['8k']=0.6690391275356948  
entropyThresholds['16k']=0.0493520022283
entropyThresholds['23k']=0.

dataFile='/Users/alomana/Google Drive File Stream/My Drive/education/ISC/2018/shared/ISB Summer Course 2018/Sandbox for afternoon leaders/case.melanoma/original/GSE72056_melanoma_single_cell_revised_v2.txt'

malignantCellsFile='/Users/alomana/scratch/malignant.%sgenes.data.csv'%resolutionLevel
malignantCellsMetadataFile='/Users/alomana/scratch/malignant.%sgenes.tumorMetadata.csv'%resolutionLevel

testingNumCells=5000 # there is a total of 4645 cells, so if this number is larger than that, all cells will be included
testingNumVar=25000 # there is a total of 23686 genes, so if this number is larger than that, all genes will be included

selectedTumors=['79','88','78','81','80','89']
selectedTumorsRanks={}
for patient in selectedTumors:
    selectedTumorsRanks[patient]=0

# 0.2. defining some variables
allGeneNames=[]
entropies=[]

# 0.3. removing previous data
if os.path.exists(malignantCellsFile) == True:
    os.remove(malignantCellsFile)

# 0.4. starting metadata files
gMetaTumor=open(malignantCellsMetadataFile,'w')
gMetaTumor.write('cell.instance,tumor.label\n')

# 1. slitting expression
with open(dataFile,'r') as f:

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

    ### dealing with trimmed versions of the data
    malignantLabels=malignantLabels[:testingNumCells]
    tumorLabels=tumorLabels[:testingNumCells]
    ####
    
    gT=open(malignantCellsFile,'a')

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
            selectedTumorsRanks[tumorLabels[i]]=selectedTumorsRanks[tumorLabels[i]]+1

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
            if malignantLabels[i] == '2' and tumorLabels[i] in selectedTumors:
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
            gT.write(geneName)
            for i in range(len(expression)):
                if malignantLabels[i] == '2' and tumorLabels[i] in selectedTumors: 
                    gT.write(',')
                    gT.write(expression[i])

            gT.write('\n')
            lineCount=lineCount+1

gT.close()

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

print(selectedTumorsRanks)
