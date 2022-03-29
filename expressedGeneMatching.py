"""
Matching between differentially expressed genes retrieved from DESeq2
Marc Meynadier
"""

#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import sys
import glob
import pandas as pd
from functools import reduce


#------------------------------------------------------------------------------#
#                              Files management                                #
#------------------------------------------------------------------------------#


def getFilenames(experiment):
    os.chdir('../data/net/6_deseq2/adult') # Changing working directory to DESeq2 results
    path=os.getcwd()
    filenames = glob.glob(path + "/*"+experiment+"*.csv")
    return filenames

def listOfFiles(filenames,experiment):
    filesNamesClean=[]
    for i in filenames:
         newName=i.rsplit('/',1)[1]
         newName=newName.rsplit(experiment+'_',1)[1]
         newName=newName.rsplit('.csv',1)[0]
         filesNamesClean.append(newName)
    return filesNamesClean

def filenamesToDataframe(filenames):
    threshold_pvalue=0.05
    dfs = [pd.read_csv(filename) for filename in filenames]
    for i in range(len(dfs)):
        dfs[i].rename(columns={ dfs[i].columns[0]: "gene" }, inplace = True)
        dfs[i]=dfs[i].drop(dfs[i][dfs[i].padj>threshold_pvalue].index)
        dfs[i]=dfs[i][dfs[i]['padj'].notna()]
    return dfs

def experimentChoice():
    print("Select your type of organisms :\n\n1 : Adult\n2 : Juvenile")
    typeOrg = int(input())
    experiment = ""
    if typeOrg == 1:
        print("Select your type of experiment :\n\n1 : Preliminary samples")
        print("2 : Spatial comparison\n3 : Temporal comparison \n4 : True transplant\n5 : Garden short")
        expType = int(input())
        if expType == 1:
            experiment = "preliminarySamples"
        elif expType == 2:
            experiment = "spatialComparison"
        elif expType == 3:
            experiment = "temporalComparison"
        elif expType == 4:
            experiment = "trueTransplant"
        elif expType == 5:
            experiment = "gardenShort"
    elif typeOrg == 2:
       print("Select your type of experiment :\n\n1 : Replica 1") 
    return experiment

#------------------------------------------------------------------------------#
#                         Shared genes computation                             #
#------------------------------------------------------------------------------#


def compareGenes(filenames,experiment):
    filesNamesClean1 = listOfFiles(filenames,experiment)
    print("\nList of output files from DESeq2 \n")
    for i in range(len(filesNamesClean1)):
        filesNamesClean1[i] = str(i+1) + " : " + filesNamesClean1[i]
        print(filesNamesClean1[i])
    print("\nWhich files do you want to compare ? (select two)\n")
    file1=int(input()) ; file2=int(input())
    df = filenamesToDataframe(filenames) 
    geneNames=list(reduce(set.intersection, map(set, [df[file1-1].gene, df[file2-1].gene])))
    lfcValuesFile1 = [] ; lfcValuesFile2 = []
    padjValuesFile1 = [] ; padjValuesFile2 = []
    dfFile1 = df[file1-1] ; dfFile2 = df[file2-1]
    for i in range(len(geneNames)):
        lfcValuesFile1.append(dfFile1['log2FoldChange'][dfFile1['gene']==geneNames[i]].values[0])
        lfcValuesFile2.append(dfFile2['log2FoldChange'][dfFile2['gene']==geneNames[i]].values[0])
        padjValuesFile1.append(dfFile1['padj'][dfFile1['gene']==geneNames[i]].values[0])
        padjValuesFile2.append(dfFile2['padj'][dfFile2['gene']==geneNames[i]].values[0])
    for i in range(len(geneNames)):
        geneNames[i]=geneNames[i].replace('TRINITY_','')
    filesNamesClean2 = listOfFiles(filenames,experiment)
    dic = {'genes':geneNames,'lfc_'+filesNamesClean2[file1-1]:lfcValuesFile1,'lfc_'+filesNamesClean2[file2-1]:lfcValuesFile2,
    'p-adj_'+filesNamesClean2[file1-1]:padjValuesFile1,'p-adj_'+filesNamesClean2[file2-1]:padjValuesFile2}
    outputDf = pd.DataFrame(dic)
    outputDf = outputDf.sort_values(by='lfc_'+filesNamesClean2[file1-1],ascending=False)
    outputDf = outputDf.reset_index(drop=True)
    print(outputDf)
    outputDf.to_csv(filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_comparison.csv',encoding='utf-8')
    

#------------------------------------------------------------------------------#
#                              Menu functions                                  #
#------------------------------------------------------------------------------#


def menu_display():
    print("\n")
    print("--------------------------------------------")
    print("|                                          |")
    print("|          Genes shared : 1                |")
    print("|                                          |")
    print("|                  Exit : 2                |")
    print("|                                          |")
    print("--------------------------------------------")
    print("\n")
    return

def menu_app():
    experiment = experimentChoice()
    filenames = getFilenames(experiment)
    while True:
        menu_display()
        answer = int(input())
        if answer==1:
            compareGenes(filenames,experiment)
        elif answer==2:
            sys.exit(0)


#------------------------------------------------------------------------------#
#                                    Main                                      #
#------------------------------------------------------------------------------#


menu_app()
