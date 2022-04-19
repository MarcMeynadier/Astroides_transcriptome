"""
Comparison of overlaping genes betwen Ontologizer output files
Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import sys
from regex import F
import pandas as pd
from functools import reduce


def filenamesToDataframe(filenames):
    print("Define your threshold value (usually 0.05)")
    threshold_pvalue=float(input())
    dfs = [pd.read_csv(filename,error_bad_lines=False,sep='\t') for filename in filenames]
    for i in range(len(dfs)):
        dfs[i]=dfs[i].drop(dfs[i][dfs[i].p>threshold_pvalue].index)
    return dfs


def getFilesNames():
    os.chdir('../data/net/7_functionnalAnnotation/ontologizer/outputResults/codingDEG') # Changing working directory to DESeq2 results
    path=os.getcwd()
    filesNamesClean=[]
    filesNames = []
    for file in os.listdir(path):
        if file.endswith(".txt"):
            filesNames.append(file)
            fileName = file.split("results_",1)[1] ; fileName = fileName.split("-Parent",1)[0]
            filesNamesClean.append(fileName)
    dfs = filenamesToDataframe(filesNames)
    print("\nList of output files from Ontologizer\n")
    for i in range(len(filesNamesClean)):
        filesNamesClean[i] = str(i+1) + " : " + filesNamesClean[i]
        print(filesNamesClean[i]) 
    print("\nHow many files do you want to compare ?\n")  
    nbrFiles = int(input())
    comparedFiles = []
    print("\nWhich files do you want to compare ?\n")
    for i in range(nbrFiles):
        comparedFiles.append(dfs[int(input())-1])
    outputDf = comparedFiles[0] 
    for i in range(len(comparedFiles)):
        try:
            outputDf = outputDf.merge(comparedFiles[i+1],how='inner',on='ID')
        except IndexError:
            pass
    print(outputDf)

getFilesNames()