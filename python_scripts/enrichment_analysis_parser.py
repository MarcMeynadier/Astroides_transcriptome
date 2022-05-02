"""
Preparing files for Ontologizer input 
Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import wget
import shutil
import pandas as pd
from urllib.request import urlopen


#------------------------------------------------------------------------------#
#                              Files management                                #
#------------------------------------------------------------------------------#

def getOntologyFileOntologizer():
    url = 'http://purl.obolibrary.org/obo/go.obo'
    wget.download(url)
    shutil.move('go.obo', '../../data/net/8_functionnalAnnotation/ontologizer')

def getOntologyFileGOMWU():
    url = 'http://purl.obolibrary.org/obo/go.obo'
    wget.download(url)
    shutil.move('go.obo', '../../data/net/8_functionnalAnnotation/GO_MWU')


def swap_columns(df, col1, col2):
    col_list = list(df.columns)
    x, y = col_list.index(col1), col_list.index(col2)
    col_list[y], col_list[x] = col_list[x], col_list[y]
    df = df[col_list]
    return df


def getAnnotationFile():
    curedFile = []
    with open('../../data/net/8_functionnalAnnotation/hmmsearchOutput.out') as f:
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
        pfamAnnot.append(splitAnnot) ; pfamCode.append(splitCode)
    dic={'genes':geneNames,'pfam_annotation':pfamAnnot,'pfam_code':pfamCode}
    annotDf=pd.DataFrame(dic)
    annotDf = annotDf.replace(to_replace ='(_i).*', value = '', regex = True)
    return annotDf

def getAssociationFile(): 
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
    finalGoCode2 = []
    for i in finalGoCode:
        i = ",".join(map(str,i))
        finalGoCode2.append(i)
    dic = {'pfam_code':pfamCodeList,'GO_annotation':finalGoAnnot,'GO_code':finalGoCode2} 
    annotSequenceDf = getAnnotationFile()
    mergeDf = pd.DataFrame(dic)
    mergeDf = mergeDf.merge(annotSequenceDf,how='inner')
    mergeDf.drop(['GO_annotation', 'pfam_annotation','pfam_code'], axis=1, inplace=True)
    mergeDf = swap_columns(mergeDf,'GO_code','genes') 
    pathOntologizer = '../../data/net/8_functionnalAnnotation/ontologizer/'
    mergeDf.to_csv(pathOntologizer+'associationFile.ids', sep='\t', header=False, index=False)

def getPopulationFile():
    curedFile = []
    with open('../../data/net/8_functionnalAnnotation/ontologizer/adult_trinity_longest_CD_HIT.Trinity.fasta') as f:
        contents = f.readlines()
    for i in contents:
        if 'TRINITY' in i:
            i = i.split("TRINITY_",1)[1]
            i = i.split("_i",1)[0]
            curedFile.append(i)
    f=open('populationFile.txt', 'a')
    f.writelines("%s\n" % i for i in curedFile)
    f.close()
    shutil.move('populationFile.txt', '../../data/net/8_functionnalAnnotation/ontologizer')

def getProteinSequences():
    with open('../../data/net/8_functionnalAnnotation/transdecoderOutput.pep') as f:
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
        geneNames[i]=geneNames[i].replace('>','')
        geneNames[i]=geneNames[i].split(' ',1)[0]
        geneNames[i]=geneNames[i].split('_i',1)[0]  
    for i in range(len(proteinSequences)):
        proteinSequences[i] = proteinSequences[i].replace('\n','')
    dic = {'genes':geneNames,'protein_sequence':proteinSequences}
    sequencesDf = pd.DataFrame(dic)
    return sequencesDf

def getStudysetFileOntologizer(threshold):
    protDf = getProteinSequences() 
    folderOrg = 'adult'
    os.chdir('../../data/net/7_deseq2/adultTranscriptome/'+folderOrg) # Changing working directory to DESeq2 results
    path=os.getcwd()
    for file in os.listdir(path):
        if file.endswith(".csv"):
            curedFile = []
            txtName=file.replace(".csv",".txt")
            csvFile = pd.read_csv(file)
            csvFile.rename(columns={ csvFile.columns[0]: "genes" }, inplace = True)
            csvFile = csvFile[csvFile.padj<threshold] 
            csvFile = csvFile[csvFile['padj'].notna()]  
            csvFile = csvFile.merge(protDf,how='inner',on='genes') 
            geneNames = csvFile.iloc[:,0].tolist()
            for i in geneNames:
                if 'TRINITY' in i:
                    i = i.split("TRINITY_",1)[1]
                    curedFile.append(i)
            f=open(txtName,'a')
            f.writelines("%s\n" % i for i in curedFile)
            f.close()
            shutil.move(txtName, '../../../8_functionnalAnnotation/ontologizer/studySamples') 
    folderOrg = 'juvenile'
    os.chdir('../'+folderOrg) # Changing working directory to DESeq2 results
    path=os.getcwd()
    for file in os.listdir(path):
        if file.endswith(".csv"):
            curedFile = []
            txtName=file.replace(".csv",".txt")
            csvFile = pd.read_csv(file)
            csvFile.rename(columns={ csvFile.columns[0]: "genes" }, inplace = True)
            csvFile = csvFile[csvFile.padj<threshold] 
            csvFile = csvFile[csvFile['padj'].notna()]  
            csvFile = csvFile.merge(protDf,how='inner',on='genes') 
            geneNames = csvFile.iloc[:,0].tolist()
            for i in geneNames:
                if 'TRINITY' in i:
                    i = i.split("TRINITY_",1)[1]
                    curedFile.append(i)
            f=open(txtName,'a')
            f.writelines("%s\n" % i for i in curedFile)
            f.close()
            shutil.move(txtName, '../../../8_functionnalAnnotation/ontologizer/studySamples')

def getStudysetFileGOMWU(threshold):
    folderOrg = 'adult'
    os.chdir('../../data/net/7_deseq2/adultTranscriptome/'+folderOrg) # Changing working directory to DESeq2 results
    path=os.getcwd()
    pathGOMWU='../../../8_functionnalAnnotation/GO_MWU/'
    for file in os.listdir(path):
        if file.endswith(".csv"):
            curedFile = []
            fileName = file.replace('.csv','')
            csvFile = pd.read_csv(file)
            csvFile = csvFile[csvFile.padj<threshold]
            csvFile = csvFile[csvFile['padj'].notna()]
            geneNames = csvFile.iloc[:,0].tolist()
            lfcValues = csvFile.iloc[:,2].tolist() 
            for i in geneNames:
                if 'TRINITY' in i:
                    i = i.split("TRINITY_",1)[1]
                    curedFile.append(i)
            dic = {'genes':curedFile,'lfc':lfcValues}
            outputDf = pd.DataFrame(dic) 
            outputDf.to_csv(pathGOMWU+fileName+'_GOMWU.csv',encoding='utf-8',index=False)
    folderOrg = 'juvenile'
    os.chdir('../'+folderOrg) # Changing working directory to DESeq2 results
    path=os.getcwd()
    for file in os.listdir(path):
        if file.endswith(".csv"):
            curedFile = []
            fileName = file.replace('.csv','')
            csvFile = pd.read_csv(file)
            csvFile = csvFile[csvFile.padj<threshold]
            csvFile = csvFile[csvFile['padj'].notna()]
            geneNames = csvFile.iloc[:,0].tolist()
            lfcValues = csvFile.iloc[:,2].tolist() 
            for i in geneNames:
                if 'TRINITY' in i:
                    i = i.split("TRINITY_",1)[1]
                    curedFile.append(i)
            dic = {'genes':curedFile,'lfc':lfcValues}
            outputDf = pd.DataFrame(dic) 
            outputDf.to_csv(pathGOMWU+fileName+'_GOMWU.csv',encoding='utf-8',index=False)