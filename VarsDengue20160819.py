# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 11:36:15 2016

@author: iryna
"""

#### this script requires transcripts/features file with lines=features columns=patients. A phenotype file with 2 phenotypes coded 1 and 0. 

usrname='iryna'
pswd=


dirShared ="/Users/iryna/Documents/isotonic/SharedBennoIryna/Dengue/20160727_2009Classif/" # "/home/iryna/isotonic_regression/Dengue/20160817_2009Classif/" #


dirIsotonInputs = dirShared + "data/"
patternTranscripts=dirIsotonInputs+"*EXPRESSION*.txt" #tab delimited matrix with a header and and index: columns=patients lines=transcript expression
patternPhenotypes=dirIsotonInputs+"*PHENOTYPE*" # file where first line=#of patients, other lines are phenotypes of patients (0 or 1) put in the same order as the patients in the transcript file



fannot = "/home/iryna/isotonic_regression/Dengue/HTA-2_0_probeset_annotations.txt"

isoton="/home/iryna/isotonic_regression/20151226_isoton/20151226_isoton"

nApproxJobsRealData=1000
nApproxJobsCvTestSetData=2000
nApproxJobsCvIntSetData=3000
nApproxJobsCvIntRealData=1000




####usually you won't modify these parameters:

mathematica='math' #'/Applications/Mathematica.app/Contents/MacOS/MathKernel'
binarize='writeBinaryData_iryna_17jun16.m' #dirShared+'code/writeBinaryData_iryna_17jun16.m'


dirIsotonOutputs = dirShared + "isotonicResults/"
relativePathIntSetCv="cvInternalSets/"
relativePathTestSetCv="cvTestSet/"
relativePathReal="real/"
relativePathCvIntReal=relativePathReal+relativePathIntSetCv
SuffixCatAndSort="_concatenatedAndSorted/"
dirCvTestSet=dirIsotonInputs+relativePathTestSetCv
PrefixFCvSet="Forisoton_"
PrefixFOutCvSet="Fromisoton_"
SuffixFCvTestSet="_LOTest"
PrefixFPhenoCvSet="Phenotype_"

dirCvIntSet=dirIsotonInputs+relativePathIntSetCv
SuffixFCvIntSet="_LOInt"

dirCvIntReal=dirIsotonInputs+relativePathReal+relativePathIntSetCv


fRealScript=dirIsotonInputs+'real.sh'
fTestSetScript=dirIsotonInputs+'testSet.sh' 
fIntSetScript=dirIsotonInputs+'intSet.sh'    
fRealIntSetScript=dirIsotonInputs+'realIntSet.sh'



