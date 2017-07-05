### this pipeline brings you from FASTQ files to expression

import os,sys

def cuffnormCaller():

    '''
    this function calls cuffnorm using all generated abundance binary files by cuffquant to produce expression tables
    '''

    outputDir=cufflinksDir+'allSamples'

    term1='time cuffnorm %s '%(genomeAnnotationFile)
    term2=' '.join(abundanceFiles)+' '
    term3='-o %s '%outputDir
    term4='-p %s '%numberOfThreads
    term5='--library-type fr-firststrand '
    term6='-L '+','.join(labels)+' '
    cmd=term1+term2+term3+term4+term5+term6
    
    print 
    print cmd
    print 

    os.system(cmd)

    return None


def cuffquantCaller(inputFile):

    '''
    this function calls the first step of the cufflinks pipeline.
    '''

    label=inputFile.split('/')[-2]
    outputDir=cufflinksDir+label
    
    term1='time cuffquant %s '%(genomeAnnotationFile)
    term2='-o %s '%outputDir
    term3='-p %s '%numberOfThreads
    term4='-M %s '%maskFile
    term5='--library-type fr-firststrand '
    term6='--multi-read-correct '
    cmd=term1+term2+term3+term4+term5+term6+inputFile

    print
    print cmd
    print

    os.system(cmd)

    return None

def STARcalling(name):

    '''
    this function calls STAR
    '''
    
    finalDir=bamFilesDir+name+'/'
    if os.path.exists(finalDir) == False:
        os.mkdir(finalDir)

    read1=cleanFASTQdir+name+'.forward.clean.fastq'
    read2=cleanFASTQdir+name+'.reverse.clean.fastq'
    sampleString=read1 + ' ' + read2

    flag1=' --genomeDir %s'%genomeIndexDir
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --readFilesIn %s'%sampleString   
    flag4=' --outFileNamePrefix %s'%finalDir
    flag5=' --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 5357465103'

    cmd='time '+STARexecutable+flag1+flag2+flag3+flag4+flag5
    
    print
    print cmd
    print
    os.system(cmd)
    
    return None

#####
##### MAIN
#####

print 'welcome to our transcriptomics tutorial code.'
print 

# 0. user defined variables

# 0.1. user defined paths
localTranscriptomicsTutorialPath='/proj/omics4tb/alomana/education/transcriptomicsTutorial/' ### very important to edit this line as where you placed your folder. Mine was /proj/omics4tb/alomana/education/transcriptomicsTutorial/
FASTQdir=localTranscriptomicsTutorialPath+'data/FASTQ/'
cleanFASTQdir=localTranscriptomicsTutorialPath+'data/cleanFASTQ/'

genomeIndexDir=localTranscriptomicsTutorialPath+'data/starIndex'
genomeFastaFile=localTranscriptomicsTutorialPath+'data/annotation/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'              
genomeAnnotationFile=localTranscriptomicsTutorialPath+'data/annotation/Saccharomyces_cerevisiae.R64-1-1.34.gff3'
maskFile=localTranscriptomicsTutorialPath+'data/annotation/Saccharomyces_cerevisiae.R64-1-1.34.masked.gff3'

bamFilesDir=localTranscriptomicsTutorialPath+'data/BAM/'
cufflinksDir=localTranscriptomicsTutorialPath+'data/cufflinks/'

transcriptomeSequenceFile=localTranscriptomicsTutorialPath+'data/annotation/orf_coding.fasta'
transcriptomeIndex=localTranscriptomicsTutorialPath+'data/kallistoIndex/'
kallistoDir=localTranscriptomicsTutorialPath+'data/kallisto/'

# 0.2. user defined executable paths
trimmomaticPath='/proj/omics4tb/alomana/software/Trimmomatic-0.36/'
STARexecutable='/proj/omics4tb/alomana/software/STAR-2.5.3a/bin/Linux_x86_64/STAR'
kallistoExecutable='/proj/omics4tb/alomana/software/kallisto_linux-v0.43.1/kallisto'

# 0.3. miscellanea of user defined variables
numberOfThreads=16 # maximum number of threads in osiris is 16

# 1. cleaning reads
print 'cleaning reads...'

# 1.1. locating input files
allFiles=os.listdir(FASTQdir)
allNames=[element.split('_R')[0] for element in allFiles]
names=list(set(allNames))

# 1.2. calling Trimmomatic
for name in names:
    
    inputForward='%s%s_R1_001.fastq'%(FASTQdir,name) # syntax hint: when you want to use the value of a variable in a text string, use %s
    inputReverse='%s%s_R2_001.fastq'%(FASTQdir,name)
    
    outputForwardClean='%s%s.forward.clean.fastq'%(cleanFASTQdir,name)
    outputForwardJunk='%s%s.forward.junk.fastq'%(cleanFASTQdir,name)
    outputReverseClean='%s%s.reverse.clean.fastq'%(cleanFASTQdir,name)
    outputReverseJunk='%s%s.reverse.junk.fastq'%(cleanFASTQdir,name)

    command='time java -jar %strimmomatic-0.36.jar PE -phred33 %s %s %s %s %s %s ILLUMINACLIP:%sadapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(trimmomaticPath,inputForward,inputReverse,outputForwardClean,outputForwardJunk,outputReverseClean,outputReverseJunk,trimmomaticPath) # syntax hint: this line is way too complicated. Check the printed command on the Terminal for simpler version

    print 'about to clean %s...'%(name)
    print command
    print 
    os.system(command)
    print

# 2. command to run STAR
print
print 'mapping reads...'

# 2.1. creating genome index
print '\t creating a genome index...'

flag1=' --runMode genomeGenerate'
flag2=' --runThreadN %s'%numberOfThreads
flag3=' --genomeDir %s'%genomeIndexDir
flag4=' --genomeFastaFiles %s'%genomeFastaFile
flag5=' --sjdbGTFfile %s'%genomeAnnotationFile
flag6=' --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 75 --genomeSAindexNbases 8'
flag7=' --outFileNamePrefix %s/'%genomeIndexDir

cmd='time '+STARexecutable+flag1+flag2+flag3+flag4+flag5+flag6+flag7
    
print
print cmd
print
os.system(cmd)

# 2.2. actual mapping of reads
print '\t mapping each sample to the genome...'
for name in names:
    STARcalling(name)

# 3. quantify alignments with cufflinks
print
print 'quantifying alignments with cufflinks toolkit...'
bamFiles=[bamFilesDir+name+'/Aligned.sortedByCoord.out.bam' for name in names]
abundanceFiles=[cufflinksDir+name+'/abundances.cxb' for name in names]
labels=[name for name in names]

# 3.1. calling cuffquantCaller 
print 'calling cuffquant...'
for inputFile in bamFiles:
    cuffquantCaller(inputFile)

# 3.2 calling cuffnorm
print
print 'calling cuffnorm...'
cuffnormCaller()


# 4. running kallisto
print
print 'starting kallisto pipeline...'
# 4.1. creating an index of the transcriptome
print 'creating transcriptome index...'
cmd='time %s index --index=%stranscriptomeIndex.idx %s'%(kallistoExecutable,transcriptomeIndex,transcriptomeSequenceFile)
print
print cmd
print
os.system(cmd)

# 4.2. quantifying expression
print
print 'quantifying expression...'

for name in names:
    cmd='time %s quant --index=%stranscriptomeIndex.idx -o %s%s --bias --plaintext -t %s -b 100 --rf-stranded %s%s.forward.clean.fastq %s%s.reverse.clean.fastq'%(kallistoExecutable,transcriptomeIndex,kallistoDir,name,numberOfThreads,cleanFASTQdir,name,cleanFASTQdir,name)

    print
    print cmd
    print
    os.system(cmd)

# 4. farewell
print
print '... all done!'
