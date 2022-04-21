"""
Filtering DESeq2 results with enrichment analysis results
Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import glob
import pandas as pd

def matchingFiles():
    os.chdir('../data/net/7_functionnalAnnotation/functionnalGenesAnalysis/DESeq2_analysis')
    path=os.getcwd()
    filenames = glob.glob(path + "/*.csv")
    filesNamesClean=[]
    for i in filenames:
         newName=i.rsplit('/',1)[1]
         newName=newName.rsplit('.csv',1)[0]
         filesNamesClean.append(newName)
    print(filesNamesClean)

matchingFiles()


