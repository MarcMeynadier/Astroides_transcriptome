"""
DESeq2_X_ontologizer : Retrieval of differentially expressed and annotated genes with 
statistically significant Gene Ontology terms.
Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import pandas as pd


#------------------------------------------------------------------------------#
#                              Files management                                #
#------------------------------------------------------------------------------#


def gettingFiles(path):
    """
    Description
    -----------
    Allows the user to retrieve file names and contents.

    Parameters
    ----------
    path
        str, path to the CSV files needed for the analysis.

    Returns
    -------
    filesList
        list, contains the names of the files available in the folder.
    dfList
        list, contains the information contained in the CSV files available in the folder.
    """  

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


#------------------------------------------------------------------------------#
#                            Analysis computation                              #
#------------------------------------------------------------------------------#


def matchingFiles():
    """
    Description
    -----------
    Retrieves the set of results from DESeq2_analysis and ontologyzer_analysis, and if two 
    files have the same name, then the program explores the output file from DESeq2_analysis. 
    For each line, the program checks if the term GO is available in the ontologizer_analysis output file. 
    If it is, the index of the corresponding line is stored in a list. Once the index list is complete, 
    it is used to filter the output file of DESeq2_analysis, keeping only the lines whose index has been stored. 
    This procedure is repeated for each file whose name matches.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """ 

    DESeq2FilesPath='../../data/net/8_functionalAnnotation/functionalGenesAnalysis/DESeq2_analysis'
    ontologizerFilesPath='../ontologizer_analysis'
    DESeq2Files,DESeq2Df = gettingFiles(DESeq2FilesPath)
    ontologizerFiles,ontologizerDf = gettingFiles(ontologizerFilesPath)
    for i in range(len(DESeq2Files)):
        indexes_to_keep = []
        for j in range(len(ontologizerFiles)):
            if DESeq2Files[i] == ontologizerFiles[j]:
                for index, row in DESeq2Df[i].iterrows():
                    for index2, row2 in ontologizerDf[j].iterrows():
                        if row2['ID'] in row['GO_code']:
                            indexes_to_keep.append(index)
        unique_indexes_keep = list(set(indexes_to_keep))
        DESeq2Df[i] = DESeq2Df[i].take(unique_indexes_keep)
        DESeq2Df[i] = DESeq2Df[i].reset_index(drop=True) 
        if not DESeq2Df[i].empty:
            DESeq2Df[i].to_csv('../DESeq2_X_ontologizer/'+DESeq2Files[i]+'_filtered.csv',encoding='utf-8')       
    print('\nProcess is finished\n')                   
                     