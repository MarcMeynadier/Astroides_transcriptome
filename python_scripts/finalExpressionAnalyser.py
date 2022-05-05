"""
Final Expression Analyser

Marc Meynadier
"""
import os
import pandas as pd
import numpy as np
import glob
from DESeq2_analysis import experimentChoice 
from DESeq2_analysis import listOfFiles
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def getFilenames(experiment):
    script_dir = os.path.dirname(__file__)
    path = os.path.join(script_dir, '../../data/net/8_functionnalAnnotation/functionnalGenesAnalysis/DESeq2_X_ontologizer')
    os.chdir(path)
    filenames = glob.glob(path + "/*"+experiment+"*.csv")
    return filenames   


def filenamesToDf(filenames,experiment):
    filenamesClean = listOfFiles(filenames,experiment) 
    for i in range(len(filenamesClean)):
        filenamesClean[i] = str(i+1) + " : " + filenamesClean[i]
        print(filenamesClean[i])
    dfs = [pd.read_csv(filename,error_bad_lines=False,sep=',') for filename in filenames]
    print("\nHow many files do you want to compare ?\n")  
    nbrFiles = 3
    while nbrFiles!=0:
        try:
            nbrFiles = int(input())
            if nbrFiles == 0:
                print("\nYou must chose at least one file\n")
            else:
                break
        except ValueError:
            print("\nYou must indicate an integer number corresponding to the number of files you want to compare\n")
    comparedFiles = []
    print("\nWhich files do you want to compare ?\n")
    selectedFilesNames = []
    for i in range(nbrFiles):
        while True:
            try:
                number = int(input())-1
                comparedFiles.append(dfs[number])
                selectedFilesNames.append(filenamesClean[number])
            except ValueError:
                print("\nYou must indicate an integer number corresponding to the file you are chosing\n")
                continue
            break
    conditions = []
    for i in selectedFilesNames:
        i = i.split(': ',1)[1]
        if '_single' in i: 
            conditions.append(i.split('_single',1)[0])
        elif '_shared' in i:
            conditions.append(i.split('_shared',1)[0])
        elif '_unshared' in i:
            conditions.append(i.split('_unshared',1)[0])
    return comparedFiles,conditions


def protFamilies():
    catabolic_enzyme_l=['hydrolase','hydrolysis','protease','peptidase','endopeptidase','metalloendopeptidase','proteolysis','phosphatase','phospholipase'
    ,'catalytic activity','metallopeptidase','metallocarboxypeptidase','lipase','catabolic process']
    isomerase_enzyme_l=['lyase','isomerase','transferase','mutase','oxydoreductase','dehydrogenase','kinase']
    binding_enzyme_l=['carbohydrate binding','chitin binding','ion binding','acid binding','protein binding']
    ion_regulation_l=['proton','ion transport','calcium ion','zinc ion','magnesium ion','channel','transmembrane transport','intermembrane','translocase',
    'iron ion','copper ion','cation transport','anion transport','sulfate transport']
    genetic_regulation_l=['DNA','histone','nucleic acid','telomere','helicase']
    transcriptomic_regulation_l=['translation','ribosome','RNA','elongation']
    cytoskeleton_regulation_l = ['microtubule','cytoskeleton','mitotic spindle','cytokinesis','extracellular matrix','actin']
    redox_regulation_l = ['oxydoreductase','oxygenase','antioxydant','redox homeostasis','NAD','NADH','dehydrogenase','FAD','FADH']
    immune_regulation_l = ['tumor','necrosis','apoptosis','immune']
    other_l = ['extracellular','signal transduction','vesicle','lipid transport']
    protFam = {'catabolic_enzyme':catabolic_enzyme_l,'isomerase_enzyme':isomerase_enzyme_l,'binding_enzyme':binding_enzyme_l,
    'ion_regulation':ion_regulation_l,'genetic_regulation':genetic_regulation_l,'translation_regulation':transcriptomic_regulation_l,
    'cytoskeleton_regulation':cytoskeleton_regulation_l,'redox_regulation':redox_regulation_l,'immune_regulation':immune_regulation_l,'miscellaneous_functions':other_l}
    return protFam

def protResults():
    catabolic_enzyme_l = [] ; isomerase_enzyme_l = [] ; binding_enzyme_l = [] ; ion_regulation_l = [] ; genetic_regulation_l = [] ; transcriptomic_regulation_l = []
    cytoskeleton_regulation_l = [] ; redox_regulation_l = [] ; immune_regulation_l = [] ; other_l = []
    protRes = {'catabolic_enzyme':catabolic_enzyme_l,'isomerase_enzyme':isomerase_enzyme_l,'binding_enzyme':binding_enzyme_l,
    'ion_regulation':ion_regulation_l,'genetic_regulation':genetic_regulation_l,'translation_regulation':transcriptomic_regulation_l,
    'cytoskeleton_regulation':cytoskeleton_regulation_l,'redox_regulation':redox_regulation_l,'immune_regulation':immune_regulation_l,'miscellaneous_functions':other_l}
    return protRes

def sortResults(df,conditions):
    protFam = protFamilies()
    protAnnot = protResults()
    protExpr = protResults()
    for index, row in df.iterrows():   
        for j in protFam:
            for k in protFam[j]:
                if k in row['GO_annotation']:  
                    annot = row['pfam_annotation'].split(' ',1)[1]     
                    protAnnot[j].append(annot)   
                    protExpr[j].append(row['lfc_'+conditions])   
                    break
    for i in protExpr:
        if len(protExpr[i]) != 0:
            averageExpr = np.mean(protExpr[i])
            protExpr[i] = averageExpr
    return protAnnot,protExpr 

def exploitResults(dfs,conditions):
    for i in range(len(dfs)):
        protAnnot,protExpr = sortResults(dfs[i],conditions[i])
        print('\nProteins results :\n')
        print(protAnnot)
        print('\nExpression results :\n')
        print(protExpr,'\n')
        #names = list(protExpr.keys())
        
        """
        plt.bar(range(len(protExpr)),values,tick_label=names)
        plt.show()
        """

typeOrg, experiment, org= experimentChoice()
filenames = getFilenames(experiment)
dfs,conditions = filenamesToDf(filenames,experiment)
exploitResults(dfs,conditions)