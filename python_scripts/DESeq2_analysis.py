"""
Matching between differentially expressed genes retrieved from DESeq2
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
    with open('../../../8_functionnalAnnotation/transdecoderOutput.pep') as f:
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
    sequencesDf = getProteinSequences()
    curedFile = []
    with open('../../../8_functionnalAnnotation/hmmsearchOutput.out') as f:
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
    curedList = []
    with open('../../../8_functionnalAnnotation/pantherOutput.out') as f:
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
    if typeOrg == 1:
        folderOrg = "/adult"
    elif typeOrg == 2:
        folderOrg = "/juvenile"
    os.chdir('../../data/net/7_deseq2/adultTranscriptome'+folderOrg) # Changing working directory to DESeq2 results
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

def filenamesToDataframe(filenames,threshold,flagCandidate):
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


def setThreshold():
    print("\nDefine your threshold value (usually 0.05)\n")
    while True:
        try:
            threshold_pvalue=float(input())
        except ValueError:
            print("\nYou must indicate a valid float ranging from 0 to 1\n")
            continue
        break
    return threshold_pvalue


def filterByCandidate():
    print("\nDo you want to filter genes by candidate genes ? Y or N\n")
    flagCandidate = input()
    if flagCandidate == 'Y':
        return flagCandidate
    elif flagCandidate == 'N':
        return flagCandidate
    else:
        print('\nYou did not choose a valid answer, no filtration by candidate genes will occur\n')
        flagCandidate = 'N'
        return flagCandidate



def singleFile(filenames,experiment,org,threshold,flagCandidate):
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
    pathFunctionnalAnnotation='../../../8_functionnalAnnotation/functionnalGenesAnalysis/DESeq2_analysis/' 
    if not outputDf.empty:
        if flagCandidate == 'Y':
            if org == 'adult':
                outputDf.to_csv(pathFunctionnalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file-1]+'_genesCandidates'+'_singleFile.csv',encoding='utf-8')
            elif org == 'juvenile':
                outputDf.to_csv(pathFunctionnalAnnotation+org+'_'+filesNamesClean2[file-1]+'_genesCandidates'+'_singleFile.csv',encoding='utf-8') 
        else:
            if org == 'adult':
                outputDf.to_csv(pathFunctionnalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file-1]+'_singleFile.csv',encoding='utf-8')
            elif org == 'juvenile':
                outputDf.to_csv(pathFunctionnalAnnotation+org+'_'+filesNamesClean2[file-1]+'_singleFile.csv',encoding='utf-8') 
        print(outputDf)
        return outputDf
    else:
        print("\nNo results are available\n")
    


def genesUnshared(filenames,experiment,org,threshold,flagCandidate):
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
        pathFunctionnalAnnotation='../../../8_functionnalAnnotation/functionnalGenesAnalysis/DESeq2_analysis/' 
        if flagCandidate == 'Y':
            outputDf.to_csv(pathFunctionnalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_genesCandidates'+'_unshared.csv',encoding='utf-8')
        else:
            outputDf.to_csv(pathFunctionnalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_unshared.csv',encoding='utf-8') 
    else:
        print("\nNo results are available\n")

def genesShared(filenames,experiment,org,threshold,flagCandidate):
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
    pathFunctionnalAnnotation='../../../8_functionnalAnnotation/functionnalGenesAnalysis/DESeq2_analysis/'
    if not outputDf.empty:
        if flagCandidate == 'Y':
            outputDf.to_csv(pathFunctionnalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_genesCandidates'+'_shared.csv',encoding='utf-8')
        else:
            outputDf.to_csv(pathFunctionnalAnnotation+org+'_'+experiment+'_'+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_shared.csv',encoding='utf-8') 
    else:
        print("\nNo results are available\n")

