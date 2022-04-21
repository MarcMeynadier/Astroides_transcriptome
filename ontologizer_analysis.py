"""
Comparison of overlaping GO terms between Ontologizer output files
Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import pandas as pd


def filenamesToDataframe(filenames):
    print("\nDefine your threshold value (usually 0.05)\n")
    while True:
        try:
            threshold_pvalue=float(input())
        except ValueError:
            print("\nYou must indicate a valid float ranging from 0 to 1\n")
            continue
        break
    dfs = [pd.read_csv(filename,error_bad_lines=False,sep='\t') for filename in filenames]
    for i in range(len(dfs)):
        dfs[i]=dfs[i].drop(dfs[i][dfs[i].p>threshold_pvalue].index)
    return dfs

def singleFileOntologizer():
    os.chdir('../data/net/7_functionnalAnnotation/ontologizer/outputResults/codingDEG') 
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
    nbrFiles = 1
    comparedFiles = []
    print("\nWhich file do you want to functionnaly enrich ?\n")
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
    selectedFiles = []
    for i in selectedNumbers:
            selectedFiles.append(filesNamesClean2[i])       
    if not outputDf.empty:
        outputDf = outputDf[['ID','name']]
        outputDf = outputDf.rename(columns={"name": "GO_term"})
        print(outputDf)
        pathFunctionnalAnnotation='../../../functionnalGenesAnalysis/ontologizer_analysis/'
        outputFileName = selectedFiles[0]
        outputDf.to_csv(pathFunctionnalAnnotation+outputFileName+'_single_file.csv',encoding='utf-8')
    else:
        print("\nNo GO terms are available for this file\n")
   

def genesUnsharedOntologizer():
    os.chdir('../data/net/7_functionnalAnnotation/ontologizer/outputResults/codingDEG') 
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
    nbrFiles = 3
    while nbrFiles!=0 or nbrFiles!=1:
        try:
            nbrFiles = int(input())
            if nbrFiles == 0 or nbrFiles == 1:
                print("\nYou must chose at least two files\n")
            else:
                break
        except ValueError:
            print("\nYou must indicate an integer number corresponding to the number of files you want to compare\n")
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
    substractDf = comparedFiles[0]
    for i in range(len(comparedFiles)):
        try:
            outputDf = outputDf.merge(comparedFiles[i+1],how='outer',on='ID')
            substractDf = outputDf.merge(comparedFiles[i+1],how='inner',on='ID')
        except IndexError:
            pass
    selectedFiles = []
    for i in selectedNumbers:
            selectedFiles.append(filesNamesClean2[i])       
    if not outputDf.empty and len(comparedFiles)!=1:
        outputDf = outputDf[['ID','name_y']]
        outputDf = outputDf.rename(columns={"name_y": "GO_term"})
        if nbrFiles > 2:
            substractID = substractDf.iloc[:,0]
            outputDf = outputDf[outputDf.ID.isin(substractID) == False]
        outputDf = outputDf.dropna(subset=['GO_term'])
        outputDf = outputDf.reset_index(drop=True)
        if not outputDf.empty:
            print(outputDf)
            pathFunctionnalAnnotation='../../../functionnalGenesAnalysis/ontologizer_analysis/'
            outputFileName = '_X_'.join(selectedFiles)
            outputDf.to_csv(pathFunctionnalAnnotation+outputFileName+'_unshared.csv',encoding='utf-8')
        else:
            print("\nNo GO terms are overlapping between those conditions with this threshold\n")
    else:
        print("\nNo GO terms are overlapping between those conditions with this threshold\n")



def genesSharedOntologizer():
    os.chdir('../data/net/7_functionnalAnnotation/ontologizer/outputResults/codingDEG') 
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
    nbrFiles = 3
    while nbrFiles!=0 or nbrFiles!=1:
        try:
            nbrFiles = int(input())
            if nbrFiles == 0 or nbrFiles == 1:
                print("\nYou must chose at least two files\n")
            else:
                break
        except ValueError:
            print("\nYou must indicate an integer number corresponding to the number of files you want to compare\n")
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
    selectedFiles = []
    for i in selectedNumbers:
            selectedFiles.append(filesNamesClean2[i])       
    if not outputDf.empty and len(comparedFiles)!=1:
        outputDf = outputDf[['ID','name_x']]
        outputDf = outputDf.rename(columns={"name_x": "GO_term"})
        print(outputDf)
        pathFunctionnalAnnotation='../../../functionnalGenesAnalysis/ontologizer_analysis/'
        outputFileName = '_X_'.join(selectedFiles)
        outputDf.to_csv(pathFunctionnalAnnotation+outputFileName+'_shared.csv',encoding='utf-8')
    else:
        print("\nNo GO terms are overlapping between those conditions with this threshold\n")
