"""
Filtering DESeq2 results with enrichment analysis results
Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import pandas as pd

def gettingFiles(path):
    os.chdir(path)
    path=os.getcwd()
    filesList = []
    dfList = []
    for file in os.listdir(path):
        if file.endswith(".csv"):   
            txtName=file.replace('.csv','')
            filesList.append(txtName)
            csvFile = pd.read_csv(file,index_col=False)
            dfList.append(csvFile)
    return filesList, dfList

def matchingFiles():
    DESeq2FilesPath='../../data/net/8_functionnalAnnotation/functionnalGenesAnalysis/DESeq2_analysis'
    ontologizerFilesPath='../ontologizer_analysis'
    DESeq2Files,DESeq2Df = gettingFiles(DESeq2FilesPath)
    ontologizerFiles,ontologizerDf = gettingFiles(ontologizerFilesPath)
    for i in range(len(DESeq2Files)):
        indexes_to_drop = []
        for j in range(len(ontologizerFiles)):
            if DESeq2Files[i] == ontologizerFiles[j]:
                for index, row in DESeq2Df[i].iterrows():
                    for index2, row2 in ontologizerDf[j].iterrows():
                        if row2['ID'] in row['GO_code']:
                            indexes_to_drop.append(index)
        unique_indexes_drop = list(set(indexes_to_drop))
        DESeq2Df[i] = DESeq2Df[i].take(unique_indexes_drop)
        DESeq2Df[i] = DESeq2Df[i].reset_index(drop=True) 
        if not DESeq2Df[i].empty:
            DESeq2Df[i].to_csv('../DESeq2_X_ontologizer/'+DESeq2Files[i]+'_filtered.csv',encoding='utf-8')       
    print('\nProcess is finished\n')                   
                     