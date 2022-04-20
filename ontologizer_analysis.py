"""
Comparison of overlaping genes betwen Ontologizer output files
Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import pandas as pd
import pathlib


def filenamesToDataframe(filenames):
    print("\nDefine your threshold value (usually 0.05)\n")
    while True:
        try:
            threshold_pvalue=float(input())
        except ValueError:
            print("\nYou must indicate a valid float ranging from 0 to 1\n")
            continue
        break
    dfs = [pd.read_csv(filename,on_bad_lines='skip',sep='\t') for filename in filenames]
    for i in range(len(dfs)):
        dfs[i]=dfs[i].drop(dfs[i][dfs[i].p>threshold_pvalue].index)
    return dfs


def ontologizerOutputAnalysis():
    os.chdir('../data/net/7_functionnalAnnotation/ontologizer/outputResults/codingDEG') # Changing working directory to DESeq2 results
    path=os.getcwd()
    filesNamesClean1=[]
    filesNames = []
    for file in os.listdir(path):
        if file.endswith(".txt"):
            filesNames.append(file)
            fileName = file.split("results_",1)[1] ; fileName = fileName.split("-Parent",1)[0]
            filesNamesClean1.append(fileName)
    dfs = filenamesToDataframe(filesNames)
    print("\nList of output files from Ontologizer\n")
    filesNamesClean2 = filesNamesClean1.copy()
    for i in range(len(filesNamesClean1)):
        filesNamesClean1[i] = str(i+1) + " : " + filesNamesClean1[i]
        print(filesNamesClean1[i]) 
    print("\nHow many files do you want to compare ?\n")  
    while True:
        try:
            nbrFiles = int(input())
        except ValueError:
            print("\nYou must indicate an integer number corresponding to the number of files you want to compare\n")
            continue
        break
    if nbrFiles == 0:
        print("You must chose at least one file")
        nbrFiles = int(input)
    comparedFiles = []
    print("\nWhich files do you want to compare ?\n")
    selectedNumbers = []
    for i in range(nbrFiles):
        while True:
            try:
                number = int(input())-1
                comparedFiles.append(dfs[number])
                selectedNumbers.append(number)
            except ValueError:
                print("\nYou must indicate an integer number corresponding to the file you are chosing\n")
                continue
            break
    outputDf = comparedFiles[0] 
    for i in range(len(comparedFiles)):
        try:
            outputDf = outputDf.merge(comparedFiles[i+1],how='inner',on='ID')
        except IndexError:
            pass
    if not outputDf.empty:
        outputDf = outputDf[['ID','name_x']]
        outputDf = outputDf.rename(columns={"name_x": "GO_term"})
        print(outputDf)
        pathFunctionnalAnnotation='../../../functionnalGenesAnalysis/ontologizer_analysis/'
        selectedFiles = []
        for i in selectedNumbers:
            selectedFiles.append(filesNamesClean2[i])
        outputFileName = '_X_'.join(selectedFiles)
        outputDf.to_csv(pathFunctionnalAnnotation+outputFileName+'.csv',encoding='utf-8')
    else:
        print("\nNo genes are overlapping between those conditions with this threshold\n")

