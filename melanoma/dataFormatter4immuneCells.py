# this script cleans original data into two files, one for malignant, other for non-malignant cells

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
dataFile='data/original/GSE72056_melanoma_single_cell_revised_v2.txt'

immuneCellsFile='data/formatted/immuneCells.xxk.genes.data.csv'
immuneCellsMetadataFile='data/formatted/immuneCells.xxk.genes.metadata.csv'
tumorCellsMetadataFile='data/formatted/tumorCells.xxk.genes.metadata.csv'

testingNumCells=5000 # there is a total of 4645 cells, so if this number is larger than that, all cells will be included
testingNumVar=25000 # there is a total of 23686 genes, so if this number is larger than that, all genes will be included

entropyThreshold=1.17875354871   # 2k genes
entropyThreshold=0.858869487712  # 4k genes
entropyThreshold=0.467954173009  # 8k genes
entropyThreshold=0.0549282251165 # 16k genes
entropyThreshold=0.              # 23k genes 

# 0.2. defining some variables
allGeneNames=[]
entropies=[]

# 0.3. removing previous data
if os.path.exists(immuneCellsFile) == True:
    os.remove(immuneCellsFile)

if os.path.exists(immuneCellsMetadataFile) == True:
    os.remove(immuneCellsMetadataFile)

if os.path.exists(tumorCellsMetadataFile) == True:
    os.remove(tumorCellsMetadataFile)

# 0.4. starting metadata files
gMetaTumor=open(tumorCellsMetadataFile,'w')
gMetaTumor.write('cell.instance,tumor.label\n')

gMetaImmune=open(immuneCellsMetadataFile,'w')
gMetaImmune.write('cell.instance,immune.label\n')

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
        if malignantLabels[i] == '1':
            gI.write(',')
            gI.write(cellIDs[i])
    gI.write('\n')

    # adding tumor cell labels
    for i in range(len(cellIDs)):
        if malignantLabels[i] == '1': 
            string2Write=cellIDs[i]+','+'Mel'+tumorLabels[i]+'\n'
            gMetaTumor.write(string2Write)
    
    # adding immune cell type
    for i in range(len(cellIDs)):
        if malignantLabels[i] == '1':
            string2Write=cellIDs[i]+','+immuneLabels[i]+'\n'
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
            if malignantLabels[i] == '1':
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
        expression=expression[:testingNumCells]
        if lineCount == testingNumVar:
            break
        ###
        
        # writing the gene name
        if duplicate == False and s > entropyThreshold:
            gI.write(geneName)
            for i in range(len(expression)):
                if malignantLabels[i] == '1': 
                    gI.write(',')
                    gI.write(expression[i])
            gI.write('\n')
            lineCount=lineCount+1

gI.close()

# 2. some prints about entropy thresholds
print('printing entropies...')
entropies.sort(reverse=True)

size=2000-1
threshold=entropies[size] # 1.17875354871
print(size,threshold)

size=4000-1
threshold=entropies[size] # 0.858869487712
print(size,threshold)

size=8000
threshold=entropies[size] # 0.467954173009
print(size,threshold)

size=16000
threshold=entropies[size] # 0.0548891548615
print(size,threshold)


