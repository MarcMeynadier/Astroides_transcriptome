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
from urllib.request import urlopen
from functools import reduce


#------------------------------------------------------------------------------#
#                              Files management                                #
#------------------------------------------------------------------------------#


def getProteinSequences():
    with open('../../../7_functionnalAnnotation/transdecoderOutput.pep') as f:
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
        geneNames[i]=geneNames[i].replace('>TRINITY_','')
        geneNames[i]=geneNames[i].split(' ',1)[0] 
    for i in range(len(proteinSequences)):
        proteinSequences[i] = proteinSequences[i].replace('\n','')
    dic = {'genes':geneNames,'protein_sequence':proteinSequences}
    sequencesDf = pd.DataFrame(dic)
    return sequencesDf

def getAnnotationFile():
    sequencesDf = getProteinSequences()
    curedFile = []
    with open('../../../7_functionnalAnnotation/hmmsearchOutput.out') as f:
        contents = f.readlines()
    for i in contents:
        if 'TRINITY' in i:
            curedFile.append(i)
    geneNames=[] ; pfamAnnot= [] ; pfamCode = []
    separator = "PF"
    for i in range(len(curedFile)):
        geneNames.append(curedFile[i].split(" - ",1)[0])
        geneNames[i]=geneNames[i].replace('TRINITY_','')
        splitDash = curedFile[i].split(" - ",1)[1]
        splitAnnot = splitDash.split(separator,1)[0] ; splitAnnot = splitAnnot.strip()
        splitCode = splitDash.split(separator,1)[1] ; splitCode = splitCode.split(".",1)[0] ; splitCode = separator + splitCode
        #splitSpace = splitPF.split() ; print(splitSpace)
        pfamAnnot.append(splitAnnot) ; pfamCode.append(splitCode)
    dic={'genes':geneNames,'pfam_annotation':pfamAnnot,'pfam_code':pfamCode}
    mergeDf=pd.DataFrame(dic)
    mergeDf = mergeDf.merge(sequencesDf,how='left')
    mergeDf = mergeDf.replace(to_replace ='(_i).*', value = '', regex = True)
    #mergeDf = mergeDf.assign(count=(mergeDf["protein_sequence"].str.len())).groupby('genes').max().drop('count',axis=1) 
    #print(mergeDf)
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
                #pass
        except IndexError:
            pass           
    finalGoAnnot.append(gatherGoAnnot[-1]) ; finalGoCode.append(gatherGoCode[-1]) 
    dic = {'pfam_code':pfamCodeList,'GO_annotation':finalGoAnnot,'GO_code':finalGoCode} 
    annotSequenceDf = getAnnotationFile()
    mergeDf = pd.DataFrame(dic)
    mergeDf = mergeDf.merge(annotSequenceDf,how='inner')
    return mergeDf

def getFilenames(typeOrg,experiment):
    if typeOrg == 1:
        folderOrg = "/adult"
    elif typeOrg == 2:
        folderOrg = "/juvenile"
    os.chdir('../data/net/6_deseq2/larvaeJuvenileAdultTranscriptome'+folderOrg) # Changing working directory to DESeq2 results
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
            experiment = "juvenile" 
    return typeOrg, experiment


#------------------------------------------------------------------------------#
#                         Shared genes computation                             #
#------------------------------------------------------------------------------#


def singleFile(filenames,experiment):
    filesNamesClean = listOfFiles(filenames,experiment)
    for i in range(len(filesNamesClean)):
        filesNamesClean[i] = str(i+1) + " : " + filesNamesClean[i]
        print(filesNamesClean[i])
    print("\nWhich file do you want to annotate ?\n")
    file=int(input())
    df = filenamesToDataframe(filenames) 
    geneNames = list(df[file-1].gene) 
    lfcValuesFile = [] 
    padjValuesFile = [] 
    dfFile = df[file-1] 
    for i in range(len(geneNames)):
        lfcValuesFile.append(dfFile['log2FoldChange'][dfFile['gene']==geneNames[i]].values[0])
        padjValuesFile.append(dfFile['padj'][dfFile['gene']==geneNames[i]].values[0])
    for i in range(len(geneNames)):
        geneNames[i]=geneNames[i].replace('TRINITY_','')
    filesNamesClean2 = listOfFiles(filenames,experiment)
    dic = {'genes':geneNames,'lfc_'+filesNamesClean2[file-1]:lfcValuesFile,
    'p-adj_'+filesNamesClean2[file-1]:padjValuesFile}
    outputDf = pd.DataFrame(dic)
    sequenceAnnotDf = pfam2goFile()
    outputDf = outputDf.merge(sequenceAnnotDf,how='left',on='genes')
    outputDf = outputDf.sort_values(by='lfc_'+filesNamesClean2[file-1],ascending=False)
    outputDf = outputDf.dropna(subset=['pfam_code'])
    outputDf = outputDf.drop_duplicates(subset=['pfam_code'])
    outputDf = outputDf.reset_index(drop=True) 
    print(outputDf)
    pathFunctionnalAnnotation='../../../7_functionnalAnnotation/expressedGeneMatching/' 
    outputDf.to_csv(pathFunctionnalAnnotation+filesNamesClean2[file-1]+'_'+experiment+'_single_file_annotation.csv',encoding='utf-8')


def genesUnshared(filenames,experiment):
    filesNamesClean = listOfFiles(filenames,experiment)
    for i in range(len(filesNamesClean)):
        filesNamesClean[i] = str(i+1) + " : " + filesNamesClean[i]
        print(filesNamesClean[i])
    print("\nWhich file do you want to compare ?\n(First choice : expressed genes)")
    file1=int(input()) ; file2=int(input())
    df = filenamesToDataframe(filenames) 
    geneNames = list(df[file1-1].gene) ; geneNames2 = list(df[file2-1].gene)
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
        lfcValuesFile.append(dfFile1['log2FoldChange'][dfFile1['gene']==geneNames[i]].values[0])
        padjValuesFile.append(dfFile1['padj'][dfFile1['gene']==geneNames[i]].values[0])
    for i in range(len(geneNames)):
        geneNames[i]=geneNames[i].replace('TRINITY_','')
    filesNamesClean2 = listOfFiles(filenames,experiment)
    dic = {'genes':geneNames,'lfc_'+filesNamesClean2[file1-1]:lfcValuesFile,
    'p-adj_'+filesNamesClean2[file1-1]:padjValuesFile}
    outputDf = pd.DataFrame(dic)
    sequenceAnnotDf = pfam2goFile()
    outputDf = outputDf.merge(sequenceAnnotDf,how='left',on='genes')
    outputDf = outputDf.sort_values(by='lfc_'+filesNamesClean2[file1-1],ascending=False)
    outputDf = outputDf.dropna(subset=['pfam_code'])
    outputDf = outputDf.drop_duplicates(subset=['pfam_code'])
    outputDf = outputDf.reset_index(drop=True)
    print(outputDf)
    pathFunctionnalAnnotation='../../../7_functionnalAnnotation/expressedGeneMatching/' 
    outputDf.to_csv(pathFunctionnalAnnotation+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_'+experiment+'_unshared_genes_comparison.csv',encoding='utf-8')

def genesShared(filenames,experiment):
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
    sequenceAnnotDf = pfam2goFile()
    outputDf = outputDf.merge(sequenceAnnotDf,how='left',on='genes')
    outputDf = outputDf.sort_values(by='lfc_'+filesNamesClean2[file1-1],ascending=False)
    outputDf = outputDf.dropna(subset=['pfam_code'])
    outputDf = outputDf.drop_duplicates(subset=['pfam_code'])
    outputDf = outputDf.reset_index(drop=True)
    print(outputDf)
    pathFunctionnalAnnotation='../../../7_functionnalAnnotation/expressedGeneMatching/'
    outputDf.to_csv(pathFunctionnalAnnotation+filesNamesClean2[file1-1]+"_X_"+filesNamesClean2[file2-1]+'_'+experiment+'_shared_genes_comparison.csv',encoding='utf-8')
    

#------------------------------------------------------------------------------#
#                              Menu functions                                  #
#------------------------------------------------------------------------------#


def menu_display():
    print("\n")
    print("--------------------------------------------")
    print("|                                          |")
    print("|    Single file annotation : 1            |")
    print("|                                          |")
    print("|            Genes unshared : 2            |")
    print("|                                          |")
    print("|              Genes shared : 3            |")
    print("|                                          |")
    print("|                      Exit : 4            |")
    print("|                                          |")
    print("--------------------------------------------")
    print("\n")
    return

def menu_app():
    typeOrg, experiment = experimentChoice()
    filenames = getFilenames(typeOrg,experiment)
    while True:
        menu_display()
        answer = int(input())
        if answer==1:
            singleFile(filenames,experiment)
        elif answer==2:
            genesUnshared(filenames,experiment)
        elif answer==3:
            genesShared(filenames,experiment)
        elif answer==4:
            sys.exit(0)


#------------------------------------------------------------------------------#
#                                    Main                                      #
#------------------------------------------------------------------------------#


menu_app()

     