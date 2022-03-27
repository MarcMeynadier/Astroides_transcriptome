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
    output=list(reduce(set.intersection, map(set, [df[file1-1].gene, df[file2-1].gene])))
    for i in range(len(output)):
        output[i]=output[i].replace('TRINITY_','')
    outputStr =", ".join(str(elem) for elem in output)
    print("\nShared genes list : \n\n",outputStr) ; print("\n\nNumber of shared genes :\n\n",len(output),"\n")
    filesNamesClean2 = listOfFiles(filenames,experiment)
    with open(filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_comparison.txt', 'w') as f:
        f.write("Shared genes list : \n\n"+outputStr)
        f.write("\n\nNumber of shared genes :\n\n"+str(len(output)))


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

def menu_app(experiment):
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


menu_app("trueTransplant")
