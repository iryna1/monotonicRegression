# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:44:02 2016

@author: iryna
"""
fVars="VarsDengue20160919.py"
from fVars import *

import pandas as pd
import os
import glob  # to list files with specific names
import numpy as np
from subprocess import call
import math
import paramiko #to execute ssh
from sklearn import metrics
import scipy
#### this script requires transcripts/features file with lines=features columns=patients. A phenotype file with 2 phenotypes coded 1 and 0. 




"""Def functions"""
def matchNamePattern(pattern):
    lFiles=glob.glob(pattern)
    assert (len(lFiles)==1), "The expression pattern " + pattern+ " doesn't match exactely 1 element. It matches "+ str(len(lFiles))+" elements."

    nameFile=lFiles[0]
    return nameFile

def catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fNameForMathema):
    dirOutCat=dirToCat[:(len(dirToCat)-1)]+SuffixCatAndSort
    if not os.path.exists(dirOutCat):
        os.makedirs(dirOutCat)
    
    resFile=dirOutCat+fNameForMathema
    
    bashCommand="cat "+dirToCat+filePatternToCat+" > "+resFile+"_temp1.txt" 
    bashCommand2="sort -nk 1 -nk 2 "+resFile+"_temp1.txt"+" > " +resFile+"_temp2.txt"
    bashCommand3="awk 'NR==1{a[$3]=0;a[$4]=0;print $0;next}{if(!($3 in a || $4 in a)) {print $0; a[$3]=0;a[$4]=0}}' "+resFile+"_temp2.txt > "+resFile+"_temp3.txt"
    bashCommand4="head -1000 "+resFile+"_temp3.txt > "+resFile 
    
    os.system(bashCommand)    
    os.system(bashCommand2)    
    os.system(bashCommand3)
    os.system(bashCommand4)    
        
    os.remove(resFile+"_temp1.txt")
    os.remove(resFile+"_temp2.txt")
    os.remove(resFile+"_temp3.txt")
    return resFile


def createScriptFile(fScript, nFeatures, dirIsotonOutputs, relativePathSetCv, PrefixFCvSet, SuffixFCvSet, dirCvSet, PrefixFPhenoCvSet,nApproxJobsCvSetData):
    dirOutCv=dirIsotonOutputs+relativePathSetCv
    
    if not os.path.exists(dirIsotonOutputs):
            os.makedirs(dirIsotonOutputs)
    if not os.path.exists(dirOutCv):
            os.makedirs(dirOutCv)
    SetFiles=glob.glob(dirCvSet+PrefixFCvSet+"*"+SuffixFCvSet+".txt.bin")
    SetFiles.sort()
    
    
    SetPheno=glob.glob(dirCvSet+PrefixFPhenoCvSet+"*"+SuffixFCvSet+".txt")
    SetPheno.sort()
    
    
    nSetFiles=int(len(SetFiles))
    nPhenoSetFiles=int(len(SetPheno))
    assert (nSetFiles==nPhenoSetFiles), "Not the same number of phenotype and gene expression files in set directory. nPhenotypes="+str(nPhenoSetFiles)+"nSetFiles="+str(nSetFiles)+". Files are defined as: "+dirCvSet+PrefixFCvSet+"*"+SuffixFCvSet+".bin"
    
    nApproxJobsPerCvSetFile=math.ceil(nApproxJobsCvSetData/nSetFiles)
    
    
     # determining the parameter b defined such as 2^b is nFeaturesPerJob  in isoton script:
    bIsoton = math.ceil(math.log(nFeatures * (1 + math.sqrt(1 + 8 * nApproxJobsPerCvSetFile)) / (4 * nApproxJobsPerCvSetFile), 2))
        
    mIsoton = math.ceil(nFeatures / math.pow(2, bIsoton))
        
        #  nJobsRealData corresponds to the parameter we called "s" in isoton. It is equal to half of a square of length mIsoton with its diagonal)
    nJobsPerSetFile=int(((mIsoton*mIsoton +mIsoton)/2))    
    nJobsSetData=int(nJobsPerSetFile*nSetFiles)
        
    
    if os.path.isfile(fScript):
        os.remove(fScript)
           
    with open(fScript, 'a') as the_file:
        the_file.write("#!/bin/bash\n#$ -cwd \n#$ -S /bin/bash\n#$ -o /dev/null\n#$ -t 1-"+str(math.floor(nJobsSetData))+"\n\n")
        the_file.write("files=('"+"'\t'".join(SetFiles)+"')\n\n")
        the_file.write("phenotypes=('"+"'\t'".join(SetPheno)+"')\n\n")
        the_file.write("nJobsPerSetFile="+str(nJobsPerSetFile)+"\n")
        
    #file that we're deling with    
        the_file.write("nFile=$[($SGE_TASK_ID-1)/$nJobsPerSetFile]"+"\n")
    #subjob number for that file
        the_file.write("nSubJob=$[$SGE_TASK_ID-$nFile*$nJobsPerSetFile]"+"\n")
        the_file.write("file=${files[$nFile]}"+"\n")    
        the_file.write("'"+isoton+"' --exhaustiveCV -k 1  -q 2000 -b "+str(bIsoton)+" -t $nSubJob " + "$file ${phenotypes[$nFile]} > '"+dirOutCv+"'$SGE_TASK_ID."+"${file##*/}"+"'.'"+"$nSubJob"+"'.txt'"+" \n")    
    return
    
    
    




def runScriptFileOnSpice(usrname,pswd,fScript):

###ssh spice and submit to queue the fScript file analysis

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('spice', username = usrname , password = pswd)
    
    
    
    #With -sync -y , I'm waiting for all jobs to end before the command finishes.
    
    stdin, stdout, stderr= ssh.exec_command('source /etc/profile; qsub -sync y '+fScript)
    exit_status = stdout.channel.recv_exit_status()
    stdout_output = stdout.read().decode('utf8').rstrip('\n')
    print ("Stdout: ", stdout_output)
    stderr_output = stderr.read().decode('utf8').rstrip('\n')
    print ("Stderr: ", stderr_output)
    assert (exit_status == 0), "While executing "+fScript+", qsub generated the non-zero exit status: "+str(exit_status)
    return exit_status



def makeCvFilesLeaveOneOut(dftranscripts, dfphenotypes,dirOut, suffixOut,PrefixFExprCv,PrefixFPhenoCv):
    if not os.path.exists(dirOut):
        os.makedirs(dirOut)
    nColsdfTranscripts=len(dftranscripts.columns)
    for i in range (0, nColsdfTranscripts):
        keptCols=list(range(0,i))+list(range(i+1,nColsdfTranscripts))
        dfKeptCols = dftranscripts.ix[:,keptCols]
        dfPhenoKeptCols = dfphenotypes.iloc[keptCols]
        fOut=dirOut+PrefixFExprCv+str(i)+suffixOut+".txt"
        fPhenoOut=dirOut+PrefixFPhenoCv+str(i)+suffixOut+".txt"
        dfKeptCols.to_csv( fOut, sep="\t", float_format='%.3f')
        dfPhenoKeptCols.columns=[str(int(len(dfPhenoKeptCols)))]
        dfPhenoKeptCols.to_csv( fPhenoOut, float_format='%.0f', index=False, sep="\t")
        
        
def makeCvFilesLeaveThreeOut(dftranscripts, dfphenotypes,dirOut, suffixOut,PrefixFExprCv,PrefixFPhenoCv,lim):
    if not os.path.exists(dirOut):
        os.makedirs(dirOut)
        
    nCasesThisSet = int(dfphenotypes.sum())
    nControlsThisSet = int(len(dfphenotypes)-nCasesThisSet)
    
    nColsdfTranscripts=len(dftranscripts.columns)
    for i in range (0, lim):
        #if I need to change leave out patients, change here:
        keptCols=list(range(0,i))+list(range(i+1,i+nCasesThisSet))+list(range(i+nCasesThisSet+1,i+nCasesThisSet+lim))+list(range(i+nCasesThisSet+lim+1,nCasesThisSet+nControlsThisSet))
        dfKeptCols = dftranscripts.ix[:,keptCols]
        dfPhenoKeptCols = dfphenotypes.iloc[keptCols]
        fOut=dirOut+PrefixFExprCv+str(i)+"_"+str(i+nCasesThisSet)+"_"+str(i+nCasesThisSet+lim)+suffixOut+".txt"
        fPhenoOut=dirOut+PrefixFPhenoCv+str(i)+"_"+str(i+nCasesThisSet)+"_"+str(i+nCasesThisSet+lim)+suffixOut+".txt"
        dfKeptCols.to_csv( fOut, sep="\t", float_format='%.3f')
        dfPhenoKeptCols.columns=[str(int(len(dfPhenoKeptCols)))]
        dfPhenoKeptCols.to_csv( fPhenoOut, float_format='%.0f', index=False, sep="\t")
    
##def vars
ftranscripts = matchNamePattern(patternTranscripts)
nameftranscripts=ftranscripts.rsplit('/', 1)[1]
fphenotypes=matchNamePattern(patternPhenotypes)

dfphenotypes = pd.read_csv(fphenotypes, skiprows=1, names=["Phenotype"])
nCases = int(dfphenotypes.sum())
nControls = int(len(dfphenotypes)-nCases)

dftranscripts = pd.read_csv(ftranscripts, sep="\t", index_col=0)
nFeatures = len(dftranscripts)
nColsdfTranscripts = len(dftranscripts.columns)
realCols = ['s1','s2','t1','t2','dir1','dir2']
lim = len(dftranscripts.columns)-1 #int(min((nControls-1) / 2, nCases-1))
limRealIntCv=len(dftranscripts.columns)
assert (len(dftranscripts.columns)==len(dfphenotypes)), "we do not have"+str(len(dftranscripts.columns))==str(len(dfphenotypes))
"""
###Making files for cv test set: leaving 1 element out in turn

print("Starting to make cvTestSet files...")
#assert (nSetFiles==nPhenoSetFiles), "Not the same number of phenotype and gene expression files in set directory. nPhenotypes="+str(nPhenoSetFiles)+"nSetFiles="+str(nSetFiles)+". Files are defined as: "+dirCvSet+PrefixFCvSet+"*"+SuffixFCvSet+".bin"

assert (len(dftranscripts.columns) == nCases+nControls), "nCases + nControls doesn't correspont to the length of columns of dftranscripts. "
makeCvFilesLeaveOneOut (dftranscripts, dfphenotypes,dirIsotonInputs+relativePathTestSetCv, SuffixFCvTestSet,PrefixFCvSet,PrefixFPhenoCvSet)

print("files for cv test set generated")


call([mathematica, '-script', binarize, dirCvTestSet,PrefixFCvSet+"*"+".txt" ]) 
#
#
#



#####create real.sh file for cluster:

call([mathematica , '-script', binarize, dirIsotonInputs, nameftranscripts ]) # '/Users/benno/dev/mmatest/test.ma'])

print("End of binarizing. Generating scripts for cluster...")

if not os.path.exists(dirIsotonOutputs):
    os.makedirs(dirIsotonOutputs)
    
if not os.path.exists(dirIsotonOutputs+relativePathReal ):
    os.makedirs(dirIsotonOutputs+relativePathReal)
    
# determining the parameter b defined such as 2^b is nFeaturesPerJob  in isoton script:
bIsoton = math.ceil(math.log(nFeatures * (1 + math.sqrt(1 + 8 * nApproxJobsRealData)) / (4 * nApproxJobsRealData), 2))
    
mIsoton = math.ceil(nFeatures / math.pow(2, bIsoton))
    
#  nJobsRealData corresponds to the parameter we called "s" in isoton. It is equal to half of a square of length mIsoton with its diagonal)
nJobsRealData=(mIsoton*mIsoton +mIsoton)/2
        

if os.path.isfile(fRealScript):
    os.remove(fRealScript)

with open(fRealScript, 'a') as the_file:
    the_file.write("#!/bin/bash\n#$ -cwd \n#$ -S /bin/bash\n#$ -o /dev/null\n#$ -t 1-"+str(math.floor(nJobsRealData))+"\n\n")
    the_file.write(isoton+" --exhaustiveCV -k 1  -q 3000 -b "+str(bIsoton)+" -t $SGE_TASK_ID " + ftranscripts+".bin "+ fphenotypes+" > "+dirIsotonOutputs+relativePathReal+"$SGE_TASK_ID.realResults.txt \n")    
    

    
### run  real.sh 

print("Submittig real.sh to cluster...")
####ssh spice and submit to queue the real file analysis
runScriptFileOnSpice(usrname,pswd,fRealScript)






#### create and run cvTestSet.sh file for cluster
print("Submittig cvTestSet.sh to cluster...")
createScriptFile(fTestSetScript, nFeatures, dirIsotonOutputs, relativePathTestSetCv, PrefixFCvSet, SuffixFCvTestSet, dirCvTestSet, PrefixFPhenoCvSet, nApproxJobsCvTestSetData)     
runScriptFileOnSpice(usrname,pswd,fTestSetScript)




#########Preparing real file for mathematica

print("preparing files for Mathematica...")

dirToCat=dirIsotonOutputs+relativePathReal
filePatternToCat="*.txt"

fForMathema=PrefixFOutCvSet+"all_once.txt"

fForMathemaReal=catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fForMathema)

print("real file for mathematica: ready") 


    
##########Preparing Test Set files for mathematica 

for i in range (0, nColsdfTranscripts):
    subpattern=str(i)+SuffixFCvTestSet
    filePatternToCat="*"+PrefixFCvSet+subpattern+"*.txt"  
    dirToCat=dirIsotonOutputs+relativePathTestSetCv
    
    fForMathema=PrefixFOutCvSet+subpattern+"_once.txt"
    
    catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fForMathema)

print("Test Set files for mathematica: ready") 

"""
"""
#############internal cv stuff
####Making files for internal cv (to determine k): leaving 1 element out in turn for each of the previous cvs

print("Starting to make cvIntSet files...")

patternCvTestSetInput=dirIsotonInputs+relativePathTestSetCv+PrefixFCvSet+"*"+SuffixFCvTestSet+".txt"

for j in range(0,len(glob.glob(patternCvTestSetInput))):
    fCvTestSetInput=dirIsotonInputs+relativePathTestSetCv+PrefixFCvSet+str(j)+SuffixFCvTestSet+".txt"
    
    fCvTestSetPhenoInput=dirIsotonInputs+relativePathTestSetCv+PrefixFPhenoCvSet+str(j)+SuffixFCvTestSet+".txt"
    dfCvTestSetInput=pd.read_csv(fCvTestSetInput, sep="\t", index_col=0)
    dfCvTestSetPhenoInput=pd.read_csv(fCvTestSetPhenoInput,skiprows=1, names=["Phenotype"])
    prefixExprIntCv=PrefixFCvSet+str(j)+SuffixFCvTestSet+"_"
    prefixPhenoIntCv=PrefixFPhenoCvSet+str(j)+SuffixFCvTestSet+"_"
    
    makeCvFilesLeaveOneOut (dfCvTestSetInput, dfCvTestSetPhenoInput,dirIsotonInputs+relativePathIntSetCv, SuffixFCvIntSet,prefixExprIntCv,prefixPhenoIntCv)
    print (j)


   




####Making files for int  cv real data (no left out patient for testing)

print("Starting to make cvRealIntSet files...")


if not os.path.exists(dirIsotonInputs+relativePathReal):
        os.makedirs(dirIsotonInputs+relativePathReal)
makeCvFilesLeaveOneOut (dftranscripts, dfphenotypes, dirIsotonInputs+relativePathCvIntReal, SuffixFCvIntSet,PrefixFCvSet,prefixPhenoIntCv)


print("Ended cvRealIntSet generation. Binarizing...")
call([mathematica, '-script', binarize, dirIsotonInputs+relativePathCvIntReal,PrefixFCvSet+"*"+".txt" ])
call([mathematica, '-script', binarize, dirCvIntSet,PrefixFCvSet+"*"+".txt" ])

####create and run fIntSetScript file for the cluster
print("Submitting cvIntSet.sh to cluster...")
createScriptFile(fIntSetScript, nFeatures, dirIsotonOutputs, relativePathIntSetCv, PrefixFCvSet, SuffixFCvIntSet, dirCvIntSet, PrefixFPhenoCvSet, nApproxJobsCvIntSetData)     
runScriptFileOnSpice(usrname,pswd,fIntSetScript)

####create and run fRealIntSetScript file for the cluster
print("Submitting cvIntSet.sh to cluster...")
createScriptFile(fRealIntSetScript, nFeatures, dirIsotonOutputs, relativePathReal+relativePathIntSetCv, PrefixFCvSet, SuffixFCvIntSet, dirCvIntReal, PrefixFPhenoCvSet, nApproxJobsCvIntRealData)     
runScriptFileOnSpice(usrname,pswd,fRealIntSetScript)



#####preparing Internal Set files for mathematica:

dirToCat=dirIsotonOutputs+relativePathIntSetCv


for j in range (0,nColsdfTranscripts):
    
    for i in range(0, lim):
        subPattern=str(j)+SuffixFCvTestSet+"_"+str(i)+SuffixFCvIntSet
        filePatternToCat="*"+PrefixFCvSet+subPattern+"*.txt"
        fForMathema=PrefixFOutCvSet+subPattern+"_once.txt"
        catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fForMathema)
        
        
print("Internal Set files for mathematica: ready")        
        

 


###Preparing real intSet files for mathematica


dirToCat=dirIsotonOutputs+relativePathReal+relativePathIntSetCv

      
for i in range(0, limRealIntCv):
    subPattern=str(i)+SuffixFCvIntSet
    filePatternToCat="*"+PrefixFCvSet+subPattern+"*.txt"

    patternfOut=PrefixFOutCvSet+subPattern
    
    
    fForMathema=patternfOut+"_once.txt"
    catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fForMathema)
        
print("Real Internal Set files for mathematica: ready")        
    
    



#####Annotating real file: not finished


dfOutRealOnce=pd.read_csv(fForMathemaReal, sep=" ", header=None, names=realCols) #


def annotate(dfreal, ftranscripts, fannot, outannotreal):#annotate real result:
    
    dftranscripts= pd.read_csv(ftranscripts, sep="\t", header=0, index_col=0)
    


    dfannot=pd.read_csv(fannot, sep="\t", header=0)#, index_col='Probeset_ID')
    
    dfreal['t1ProbeID']=dftranscripts.index.values[dfreal['t1']]
    dfreal['t2ProbeID']=dftranscripts.index.values[dfreal['t2']]
    
    
    colDfAnnotSaved=dfannot.columns
    dfannot.columns =['t1ProbeID','t1EntrezID','t1GeneSymbol']
    dfreal=pd.merge(dfreal,dfannot, on=['t1ProbeID'])
    dfannot.columns =['t2ProbeID','t2EntrezID','t2GeneSymbol']
    dfreal=pd.merge(dfreal,dfannot, on=['t2ProbeID'])
    dfannot.columns=colDfAnnotSaved         
    
    
    dfreal.to_csv(outannotreal, sep="\t", columns=['s1','s2','t1GeneSymbol', 't2GeneSymbol',
    't1ProbeID', 't2ProbeID','t1','t2','t1EntrezID', 't2EntrezID','dir1','dir2'])      
    return       
 
annotate( dfOutRealOnce, ftranscripts, fannot, fAnnotOutRealOnce)


"""

predicted="phenotypeProbabilities_10.tsv"

fpredicted=dirShared+predicted

predictPhenotype="calculateFinalErrorWComments_kFixed_forScript.m"
lbestK=[10]


def matchNamePattern(pattern):
    lFiles=glob.glob(pattern)
    assert (len(lFiles)==1), "The expression pattern " + pattern+ " doesn't match exactely 1 element. It matches "+ str(len(lFiles))+" elements."

    nameFile=lFiles[0]
    return nameFile


###action


for bestK in lbestK:


   # call([mathematica, '-script', str(predictPhenotype), str(dirShared),str(bestK) ]) 
    
       
    fphenotypes=matchNamePattern(patternPhenotypes)
    dfphenotypes=pd.read_csv(fphenotypes, skiprows=1, header=None)
    dfpredicted=pd.read_csv(fpredicted, header=None)
    
    
    
    roc_auc=metrics.roc_auc_score(dfphenotypes, dfpredicted)

    fpr, tpr, thresholds=metrics.roc_curve(dfphenotypes, dfpredicted)
    
    
    
    
    # Plot ROC curve
    plt.plot(fpr, tpr, label='ROC curve (area = %0.3f)' % roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')  # random predictions curve
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0,     1.0])
    plt.xlabel('False Positive Rate or (1 - Specifity)')
    plt.ylabel('True Positive Rate or (Sensitivity)')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    print("k: ", bestK) 
    print("auc =", roc_auc)
       
    cat=pd.concat([dfpredicted,dfphenotypes], axis=1, )
    cat.columns=["PredictedPhenoProba","realPheno" ]
    catsort=cat.sort_values(by="PredictedPhenoProba")
    precision, recall, thresholds=metrics.precision_recall_curve(catsort["realPheno"], catsort["PredictedPhenoProba"], pos_label=1)
    
    
    cat=pd.concat([pd.DataFrame(precision),pd.DataFrame(recall)], axis=1 )
    cat.columns=["precision","recall" ]
    catsort=cat.sort_values(by="recall")
    
    
    plt.plot( catsort["recall"],catsort["precision"],label='Precision Recall curve (area = %0.3f)' % metrics.auc(catsort["recall"],catsort["precision"]))
    #plt.plot([0, 1], [0, 1], 'k--')  # random predictions curve
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision Recall curve')
    plt.legend(loc="upper right")   
    a=metrics.auc(catsort["recall"],catsort["precision"])
    
    print("aupr =", a)
    plt.savefig('/Users/Iryna/Desktop/myimage.pdf', format="pdf", dpi=1200)
    out=pd.DataFrame([roc_auc,a], index=["auc", "aupr"], columns=[bestK])
    out.to_csv(dirShared+"PerfMeasures_k="+str(bestK)+".txt" , sep= "\t")
   
