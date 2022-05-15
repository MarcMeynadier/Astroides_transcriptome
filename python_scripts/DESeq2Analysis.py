"""
DESeq2_analyser : Treatments applied to DESeq2 contrast results.

Marc Meynadier
"""

#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import glob
import pandas as pd
from urllib.request import urlopen
from functools import reduce
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


#------------------------------------------------------------------------------#
#                              Files management                                #
#------------------------------------------------------------------------------#


def getProteinSequences():
    """
    Description
    -----------
    Retrieves the output file from TransDecoder.Predict and parse it into a dataframe containing 
    the gene IDs and their associated protein sequences.

    Parameters
    ----------
    None

    Returns
    -------
    sequencesDf
        Pandas dataframe, contains the genes IDs, and the proteins sequences associated with the genes.
    """

    with open('../../../8_functionalAnnotation/transdecoderOutput.pep') as f:
        contents = f.readlines()
    geneNames = []
    proteinSequences = []
    concatProt = []
    for i in contents:
        if 'TRINITY' not in i:
            concatProt.append(i)
        else:
            geneNames.append(i)
            sequenceProt = ''
            for j in concatProt:
                sequenceProt += ''.join(j)
            proteinSequences.append(sequenceProt) 
            concatProt = []
    for i in range(len(geneNames)):
        geneNames[i]=geneNames[i].replace('>TRINITY_','TRINITY_')
        geneNames[i]=geneNames[i].split(' ',1)[0] 
    for i in range(len(proteinSequences)):
        proteinSequences[i] = proteinSequences[i].replace('\n','')
    dic = {'genes':geneNames,'protein_sequence':proteinSequences}
    sequencesDf = pd.DataFrame(dic)
    return sequencesDf

def getAnnotationFile():
    """
    Description
    -----------
    Retrieves the output file from hmmsearch, parse it to retrieve the Pfam annotation and the Pfam code of the genes. 
    Retrieves the dataframe from getProteinSequences(), then the output file from hmmsearch in order to parse it to get 
    the Pfam annotation as well as the Pfam code of the genes. The two data frames are then merged together and returned 
    as a common data frame. The merge is done so that even a non-coding gene ID is returned. 

    Parameters
    ----------
    None

    Returns
    -------
    mergeDf
       Pandas dataframe, contains successively the gene IDs, the Pfam annotations, the Pfam codes and the protein sequences. 
    """

    sequencesDf = getProteinSequences()
    curedFile = []
    with open('../../../8_functionalAnnotation/hmmsearchOutput.out') as f:
        contents = f.readlines()
    for i in contents:
        if 'TRINITY' in i:
            curedFile.append(i)
    geneNames=[] ; pfamAnnot= [] ; pfamCode = []
    separator = "PF"
    for i in range(len(curedFile)):
        geneNames.append(curedFile[i].split(" - ",1)[0])
        splitDash = curedFile[i].split(" - ",1)[1]
        splitAnnot = splitDash.split(separator,1)[0] ; splitAnnot = splitAnnot.strip()
        splitCode = splitDash.split(separator,1)[1] ; splitCode = splitCode.split(".",1)[0] ; splitCode = separator + splitCode
        pfamAnnot.append(splitAnnot) ; pfamCode.append(splitCode)
    dic={'genes':geneNames,'pfam_annotation':pfamAnnot,'pfam_code':pfamCode}
    mergeDf=pd.DataFrame(dic) 
    mergeDf = mergeDf.merge(sequencesDf,how='left')
    mergeDf = mergeDf.replace(to_replace ='(_i).*', value = '', regex = True)
    return mergeDf

def pfam2goFile():
    """
    Description
    -----------
    Retrieve the pfam2go database from the Gene Ontology site, then parse this database to retrieve the GO terms 
    and GO codes associated with each Pfam code. Since a Pfam code can contain several GO terms, several successive 
    parsings are necessary to produce a dataframe containing the Pfam annotations, the GO terms and the GO codes. 
    The dataframe from getAnnotationFile() is then retrieved and merged with the previously parsed dataframe, 
    using the common Pfam code between the two dataframes.

    Parameters
    ----------
    None

    Returns
    -------
    mergeDf
       Pandas dataframe, contains successively the gene IDs, the Pfam annotations, the Pfam codes, the GO terms, the GO codes and the protein sequences. 
    """    
     
    url = 'http://current.geneontology.org/ontology/external2go/pfam2go'
    pfam2go = urlopen(url).read().decode('utf-8') 
    pfam2goList = pfam2go.split("\n")
    curedList = []
    for i in pfam2goList:
        if "Pfam:" in i:
            curedList.append(i)
    pfamCodeList = [] ; goAnnotList = [] ; goCodeList = []
    for i in curedList:
        pfamCode = i.split()[0] ; pfamCode = pfamCode.split("Pfam:",1)[1]
        goAnnot = i.split("> ",1)[1] ; goAnnot = goAnnot.split(" ;",1)[0]
        goCode = i.rsplit("; ",1)[1] 
        pfamCodeList.append(pfamCode) ; goAnnotList.append(goAnnot) ; goCodeList.append(goCode)
    gatherGoAnnot = [] ; gatherGoCode = []
    tempGoAnnot = [] ; tempGoCode = []
    for i in range(len(pfamCodeList)):
        if pfamCodeList[i] == pfamCodeList[i-1]:
            tempGoAnnot.append(goAnnotList[i-1]) ; 
            tempGoCode.append(goCodeList[i-1])
        else:
            tempGoAnnot.append(goAnnotList[i-1]) ; tempGoCode.append(goCodeList[i-1]) 
            gatherGoAnnot.append(tempGoAnnot) ; gatherGoCode.append(tempGoCode)
            tempGoAnnot = [] ; tempGoCode = []
    gatherGoAnnot.pop(0) ; gatherGoCode.pop(0)
    manualLastGoAnnot = [] ; manualLastGoAnnot.append(goAnnotList[-2]) 
    manualLastGoAnnot.append(goAnnotList[-1]) ; gatherGoAnnot.append(manualLastGoAnnot)
    manualLastGoCode = [] ; manualLastGoCode.append(goCodeList[-2])  
    manualLastGoCode.append(goCodeList[-1]) ; gatherGoCode.append(manualLastGoCode)
    tempGoAnnot = [] ; tempGoCode = []
    finalGoAnnot = [] ; finalGoCode = []
    count = 0
    for i in range(len(pfamCodeList)):
        try:
            if pfamCodeList[i] == pfamCodeList[i+1]:
                finalGoAnnot.append(gatherGoAnnot[count]) ; 
                finalGoCode.append(gatherGoCode[count])
            else:
                finalGoAnnot.append(gatherGoAnnot[count]) ; finalGoCode.append(gatherGoCode[count]) 
                count += 1
        except IndexError:
            pass           
    finalGoAnnot.append(gatherGoAnnot[-1]) ; finalGoCode.append(gatherGoCode[-1]) 
    dic = {'pfam_code':pfamCodeList,'GO_annotation':finalGoAnnot,'GO_code':finalGoCode} 
    annotSequenceDf = getAnnotationFile()
    mergeDf = pd.DataFrame(dic)
    mergeDf = mergeDf.merge(annotSequenceDf,how='inner')
    return mergeDf

def getPantherFiles():
    """
    Description
    -----------
    Retrieves the output file from the Panther tool, parse it to turn it into a dataframe containing gene IDs, Panther annotations and Panther codes. 
    The dataframe from pfam2goFile() is retrieved, then merges with the previously parsed dataframe using the gene IDs.

    Parameters
    ----------
    None

    Returns
    -------
    mergeDf
       Pandas dataframe, contains successively the gene IDs, the Panther annotations, the Panther codes, the Pfam codes, the GO terms, 
       the GO codes, the Pfam annotations and the protein sequences.  
    """   

    curedList = []
    with open('../../../8_functionalAnnotation/pantherOutput.out') as f:
        contents = f.readlines() 
    for i in contents:
        if 'TRINITY' in i:
            curedList.append(i)
    geneList = [] ; pantherCode = [] ; pantherAnnot = []
    for i in range(len(curedList)):
        gene = curedList[i].split('\t',1)[0]
        gene = gene.split('_i',1)[0]
        code = curedList[i].split('\t',1)[1] ; code = code.split('\t',1)[0]
        annot = curedList[i].split('\t',2)[2] ; annot = annot.split('\t',1)[0]
        geneList.append(gene)
        pantherCode.append(code)
        pantherAnnot.append(annot)
    dic={'genes':geneList,'panther_annotation':pantherAnnot,'panther_code':pantherCode}
    mergeDf=pd.DataFrame(dic)
    annotSequenceDf = pfam2goFile()
    mergeDf = mergeDf.merge(annotSequenceDf,how='inner')
    return mergeDf

def getFilenames(typeOrg,experiment):
    """
    Description
    -----------
    Retrieves the output file names from DESeq2 based on the type of organism and the associated experiment.

    Parameters
    ----------
    typeOrg
        int, 1 for adults and 2 for juveniles.       
    experiment
        str, contains the type of experiment retrieved with experimentChoice().

    Returns
    -------
    filenames
        list, contains the file names of the DESeq2 results associated with their paths.
    """   

    if typeOrg == 1:
        folderOrg = "/adult"
    elif typeOrg == 2:
        folderOrg = "/juvenile"
    os.chdir('../../data/net/7_deseq2/adultTranscriptome'+folderOrg) 
    path=os.getcwd()
    filenames = glob.glob(path + "/*"+experiment+"*.csv")
    return filenames

def listOfFiles(filenames,experiment):
    """
    Description
    -----------
    Gets the list provided by getFilenames() in order to parse it to get only the file names without their paths and formats. 

    Parameters
    ----------
    filenames
        list, contains the full names and paths of the DESeq2 output files.
    experiment
        str, contains the type of experiment retrieved with experimentChoice().

    Returns
    -------
    filesNamesClean
        list, contains the file names of the DESeq2 results without their associated paths.
    """   

    filesNamesClean=[]
    for i in filenames:
         newName=i.rsplit('/',1)[1]
         newName=newName.rsplit(experiment+'_',1)[1]
         newName=newName.rsplit('.csv',1)[0]
         filesNamesClean.append(newName)
    return filesNamesClean

def filenamesToDataframe(filenames,threshold,flagCandidate):
    """
    Description
    -----------
    Retrieves CSV files from DESeq2 by their names, transforms them into dataframe before entering them into a list. 
    The genes are then filtered according to the p-value and candidate gene thresholds, and the dataframe list is returned.

    Parameters
    ----------
    filenames
        list, contains the name of files parsed by listOfFiles().
    threshold
        int, contains the value of p-value threshold defined in the settings.
    flagCandidate
        str, contains the switch activating the filtering by candidate genes defined in the settings.

    Returns
    -------
    dfs
        list, dataframe containing the filtered DESeq2 results.
    """   

    dfs = [pd.read_csv(filename) for filename in filenames] 
    if flagCandidate == 'Y':
        candidateGenes = pd.read_csv('../../../../../Astroides_transcriptome/R_scripts/candidateGenes.csv')
    for i in range(len(dfs)):
        dfs[i].rename(columns={ dfs[i].columns[0]: "genes" }, inplace = True)
        dfs[i]=dfs[i].drop(dfs[i][dfs[i].padj>threshold].index)
        dfs[i]=dfs[i][dfs[i]['padj'].notna()]
        if flagCandidate == 'Y':
            dfs[i]=dfs[i].merge(candidateGenes,on='genes',how='inner')
    return dfs

def experimentChoice():
    """
    Description
    -----------
    Allows the user to select the type of organization and experiment from which the DESeq2 results are derived.

    Parameters
    ----------
    None

    Returns
    -------
    typeOrg
        int, 1 for adults and 2 for juveniles.       
    experiment
        str, contains the type of experiment.
    org
        str, adult or juvenile.
    """   

    print("\nSelect your type of organisms :\n\n1 : Adult\n2 : Juvenile\n")
    while True:
        try:
            typeOrg = int(input())
        except ValueError:
            print("\nYou must indicate an integer value ranging from 1 to 2\n")
            continue
        break
    experiment = ""
    if typeOrg == 1:
        org = 'adult'
        print("\nSelect your type of experiment :\n\n1 : Preliminary samples")
        print("2 : Spatial comparison\n3 : Temporal comparison \n4 : True transplant\n5 : Garden short\n")
        while True:
            try:
                expType= int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 5\n")
                continue
            break 
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
            experiment = "juvenile" 
            org = "juvenile" 
    return typeOrg, experiment, org


#------------------------------------------------------------------------------#
#                            Analysis computation                              #
#------------------------------------------------------------------------------#


def singleFile(filenames,experiment,org,threshold,flagCandidate):
    """
    Description
    -----------
    From the file names with their associated paths, parse them with listOfFiles() and retrieve the associated 
    DESeq2 data with filenamesToDataframe(). The user then chooses which file he wants to analyze. The latter is 
    retrieved as well as annotations dataframe with getPantherFiles(), and both are merged thanks to the IDs of the genes 
    in common. A CSV file containing all the information is then save and returned.

    Parameters
    ----------
    filenames
        list, contains the full names and paths of the DESeq2 output files.
    experiment
        str, contains the type of experiment.
    org
        str, adult or juvenile.
    threshold
        float, p-value threshold value.
    flagCandidate
        str, Y or N (yes or no).

    Returns
    -------
    outputDf
        Pandas dataframe, contains all the information from the user-selected DESeq2 file, 
        with filtered genes and functional annotations provided as well.
    """   

    filesNamesClean = listOfFiles(filenames,experiment)
    for i in range(len(filesNamesClean)):
        filesNamesClean[i] = str(i+1) + " : " + filesNamesClean[i]
        print(filesNamesClean[i])
    print("\nWhich file do you want to annotate ?\n")
    while True:
            try:
                file=int(input())
                if file not in range(len(filesNamesClean)+1):
                    print("\nYou must indicate an integer value corresponding to the file you want to analyse\n")
                    file=int(input())
            except ValueError:
                print("\nYou must indicate an integer value corresponding to the file you want to analyse\n")
                continue
            break 
    df = filenamesToDataframe(filenames,threshold,flagCandidate) 
    geneNames = list(df[file-1].genes) 
    lfcValuesFile = [] 
    padjValuesFile = [] 
    dfFile = df[file-1] 
    for i in range(len(geneNames)):
        lfcValuesFile.append(dfFile['log2FoldChange'][dfFile['genes']==geneNames[i]].values[0])
        padjValuesFile.append(dfFile['padj'][dfFile['genes']==geneNames[i]].values[0])
    filesNamesClean2 = listOfFiles(filenames,experiment)
    dic = {'genes':geneNames,'lfc_'+filesNamesClean2[file-1]:lfcValuesFile,
    'p-adj_'+filesNamesClean2[file-1]:padjValuesFile}
    outputDf = pd.DataFrame(dic)
    sequenceAnnotDf = getPantherFiles()
    outputDf = outputDf.merge(sequenceAnnotDf,how='inner',on='genes')
    outputDf = outputDf.sort_values(by='lfc_'+filesNamesClean2[file-1],ascending=False)
    outputDf = outputDf.dropna(subset=['pfam_code'])
    outputDf = outputDf.drop_duplicates(subset=['pfam_code'])
    outputDf = outputDf.drop_duplicates(subset=['genes']) 
    s = outputDf.pop('panther_code') ; outputDf = pd.concat([outputDf,s],1)
    s = outputDf.pop('pfam_code') ; outputDf = pd.concat([outputDf,s],1)
    s = outputDf.pop('GO_code') ; outputDf = pd.concat([outputDf,s],1)
    s = outputDf.pop('protein_sequence') ; outputDf = pd.concat([outputDf,s],1)
    outputDf = outputDf.reset_index(drop=True) 
    pathFunctionalAnnotation='../../../8_functionalAnnotation/functionalGenesAnalysis/DESeq2_analysis/' 
    if not outputDf.empty:
        if flagCandidate == 'Y':
            if org == 'adult':
                outputDf.to_csv(pathFunctionalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file-1]+'_genesCandidates'+'_singleFile.csv',encoding='utf-8')
            elif org == 'juvenile':
                outputDf.to_csv(pathFunctionalAnnotation+org+'_'+filesNamesClean2[file-1]+'_genesCandidates'+'_singleFile.csv',encoding='utf-8') 
        else:
            if org == 'adult':
                outputDf.to_csv(pathFunctionalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file-1]+'_singleFile.csv',encoding='utf-8')
            elif org == 'juvenile':
                outputDf.to_csv(pathFunctionalAnnotation+org+'_'+filesNamesClean2[file-1]+'_singleFile.csv',encoding='utf-8') 
        print(outputDf)
        return outputDf
    else:
        print("\nNo results are available\n")
    


def genesUnshared(filenames,experiment,org,threshold,flagCandidate):
    """
    Description
    -----------
    From the file names with their associated paths, parse them with listOfFiles() and retrieve the associated 
    DESeq2 data with filenamesToDataframe(). The genes that are not in common between the files and their information are recovered.
    A CSV file containing all those informations is then save and returned.

    Parameters
    ----------
    filenames
        list, contains the full names and paths of the DESeq2 output files.
    experiment
        str, contains the type of experiment.
    org
        str, adult or juvenile.
    threshold
        float, p-value threshold value.
    flagCandidate
        str, Y or N (yes or no).

    Returns
    -------
    outputDf
        Pandas dataframe, contains all the information from the user-selected DESeq2 files, 
        with filtered genes and functional annotations provided as well.
    """   

    filesNamesClean = listOfFiles(filenames,experiment)
    filesNamesClean2=filesNamesClean.copy()
    for i in range(len(filesNamesClean)):
        filesNamesClean[i] = str(i+1) + " : " + filesNamesClean[i]
        print(filesNamesClean[i])
    print("\nWhich file do you want to compare ?\n(First choice is the comparison benchmark)")
    while True:
            try:
                file1=int(input()) ; file2=int(input())
            except ValueError:
                print("\nYou must indicate an integer value corresponding to the files you want to analyse\n")
                continue
            break 
    df = filenamesToDataframe(filenames,threshold,flagCandidate) 
    geneNames = list(df[file1-1].genes) ; geneNames2 = list(df[file2-1].genes)
    for i in geneNames:
        flag = 0
        for j in geneNames2:
            if i == j:
                flag = 1
        if flag == 0:
            geneNames.remove(i)
    lfcValuesFile = [] 
    padjValuesFile = [] 
    dfFile1 = df[file1-1] 
    for i in range(len(geneNames)):
        lfcValuesFile.append(dfFile1['log2FoldChange'][dfFile1['genes']==geneNames[i]].values[0])
        padjValuesFile.append(dfFile1['padj'][dfFile1['genes']==geneNames[i]].values[0])
    dic = {'genes':geneNames,'lfc_'+filesNamesClean2[file1-1]:lfcValuesFile,
    'p-adj_'+filesNamesClean2[file1-1]:padjValuesFile}
    outputDf = pd.DataFrame(dic)
    sequenceAnnotDf = getPantherFiles()
    outputDf = outputDf.merge(sequenceAnnotDf,how='left',on='genes')
    outputDf = outputDf.sort_values(by='lfc_'+filesNamesClean2[file1-1],ascending=False)
    outputDf = outputDf.dropna(subset=['pfam_code'])
    outputDf = outputDf.drop_duplicates(subset=['pfam_code'])
    outputDf = outputDf.drop_duplicates(subset=['genes']) 
    s = outputDf.pop('panther_code') ; outputDf = pd.concat([outputDf,s],1)
    s = outputDf.pop('pfam_code') ; outputDf = pd.concat([outputDf,s],1)
    s = outputDf.pop('GO_code') ; outputDf = pd.concat([outputDf,s],1)
    s = outputDf.pop('protein_sequence') ; outputDf = pd.concat([outputDf,s],1)
    outputDf = outputDf.reset_index(drop=True)
    if not outputDf.empty:
        print(outputDf)
        pathFunctionalAnnotation='../../../8_functionalAnnotation/functionalGenesAnalysis/DESeq2_analysis/' 
        if flagCandidate == 'Y':
            outputDf.to_csv(pathFunctionalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_genesCandidates'+'_unshared.csv',encoding='utf-8')
        else:
            outputDf.to_csv(pathFunctionalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_unshared.csv',encoding='utf-8') 
    else:
        print("\nNo results are available\n")

def genesShared(filenames,experiment,org,threshold,flagCandidate):
    """
    Description
    -----------
    From the file names with their associated paths, parse them with listOfFiles() and retrieve the associated 
    DESeq2 data with filenamesToDataframe(). The genes that are in common between the files and their information are recovered.
    A CSV file containing all those informations is then save and returned.

    Parameters
    ----------
    filenames
        list, contains the full names and paths of the DESeq2 output files.
    experiment
        str, contains the type of experiment.
    org
        str, adult or juvenile.
    threshold
        float, p-value threshold value.
    flagCandidate
        str, Y or N (yes or no).

    Returns
    -------
    outputDf
        Pandas dataframe, contains all the information from the user-selected DESeq2 files, 
        with filtered genes and functional annotations provided as well.
    """   

    filesNamesClean1 = listOfFiles(filenames,experiment)
    filesNamesClean2=filesNamesClean1.copy()
    print("\nList of output files from DESeq2 \n")
    for i in range(len(filesNamesClean1)):
        filesNamesClean1[i] = str(i+1) + " : " + filesNamesClean1[i]
        print(filesNamesClean1[i])
    print("\nWhich files do you want to compare ? (select two)\n")
    while True:
            try:
                file1=int(input()) ; file2=int(input())
            except ValueError:
                print("\nYou must indicate an integer value corresponding to the files you want to analyse\n")
                continue
            break 
    df = filenamesToDataframe(filenames,threshold,flagCandidate) 
    geneNames=list(reduce(set.intersection, map(set, [df[file1-1].genes, df[file2-1].genes])))
    lfcValuesFile1 = [] ; lfcValuesFile2 = []
    padjValuesFile1 = [] ; padjValuesFile2 = []
    dfFile1 = df[file1-1] ; dfFile2 = df[file2-1]
    for i in range(len(geneNames)):
        lfcValuesFile1.append(dfFile1['log2FoldChange'][dfFile1['genes']==geneNames[i]].values[0])
        lfcValuesFile2.append(dfFile2['log2FoldChange'][dfFile2['genes']==geneNames[i]].values[0])
        padjValuesFile1.append(dfFile1['padj'][dfFile1['genes']==geneNames[i]].values[0])
        padjValuesFile2.append(dfFile2['padj'][dfFile2['genes']==geneNames[i]].values[0])
    dic = {'genes':geneNames,'lfc_'+filesNamesClean2[file1-1]:lfcValuesFile1,'lfc_'+filesNamesClean2[file2-1]:lfcValuesFile2,
    'p-adj_'+filesNamesClean2[file1-1]:padjValuesFile1,'p-adj_'+filesNamesClean2[file2-1]:padjValuesFile2}
    outputDf = pd.DataFrame(dic)
    sequenceAnnotDf = getPantherFiles()
    outputDf = outputDf.merge(sequenceAnnotDf,how='left',on='genes')
    outputDf = outputDf.sort_values(by='lfc_'+filesNamesClean2[file1-1],ascending=False)
    outputDf = outputDf.dropna(subset=['pfam_code'])
    outputDf = outputDf.drop_duplicates(subset=['pfam_code'])
    outputDf = outputDf.drop_duplicates(subset=['genes']) 
    s = outputDf.pop('panther_code') ; outputDf = pd.concat([outputDf,s],1)
    s = outputDf.pop('pfam_code') ; outputDf = pd.concat([outputDf,s],1)
    s = outputDf.pop('GO_code') ; outputDf = pd.concat([outputDf,s],1)
    s = outputDf.pop('protein_sequence') ; outputDf = pd.concat([outputDf,s],1)
    outputDf = outputDf.reset_index(drop=True)
    print(outputDf)
    pathFunctionalAnnotation='../../../8_functionalAnnotation/functionalGenesAnalysis/DESeq2_analysis/'
    if not outputDf.empty:
        if flagCandidate == 'Y':
            outputDf.to_csv(pathFunctionalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_genesCandidates'+'_shared.csv',encoding='utf-8')
        else:
            outputDf.to_csv(pathFunctionalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_shared.csv',encoding='utf-8') 
    else:
        print("\nNo results are available\n")

