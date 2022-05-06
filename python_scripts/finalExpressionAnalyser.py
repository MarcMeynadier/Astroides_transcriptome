"""
Final Expression Analyser

Marc Meynadier
"""

import os
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def getFilenamesFinal(experiment,threshold,flagCandidate):
    script_dir = os.path.dirname(__file__)
    path = os.path.join(script_dir, '../../data/net/8_functionnalAnnotation/functionnalGenesAnalysis/DESeq2_X_ontologizer')
    os.chdir(path)
    filenames = []
    files = glob.glob('*'+experiment+'*.csv')
    if flagCandidate == 'Y':
        for f in files : 
            if 'genesCandidates' in f : 
                filenames.append(f)
    else:
        for f in files : 
            if 'genesCandidates' not in f : 
                filenames.append(f)
    return filenames


def filenamesToDfFinal(filenames,experiment,flagCandidate):
    filenamesClean = []
    for i in filenames:
        if flagCandidate == 'Y':
           i = i.replace('_genesCandidates','')
        i = i.split(experiment+'_',1)[1]
        i = i.split('_filtered',1)[0]
        filenamesClean.append(i.split('.csv',1)[0]) 
    for i in range(len(filenamesClean)):
        filenamesClean[i] = str(i+1) + " : " + filenamesClean[i]
        print(filenamesClean[i])
    dfs = [pd.read_csv(filename,error_bad_lines=False,sep=',') for filename in filenames]
    if len(filenamesClean) == 0:
        return dfs,filenamesClean
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
    return comparedFiles,conditions,experiment


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
    protFam = {'catabolic enzyme':catabolic_enzyme_l,'isomerase enzyme':isomerase_enzyme_l,'binding enzyme':binding_enzyme_l,
    'ion regulation':ion_regulation_l,'genetic regulation':genetic_regulation_l,'translation regulation':transcriptomic_regulation_l,
    'cytoskeleton regulation':cytoskeleton_regulation_l,'redox regulation':redox_regulation_l,'immune regulation':immune_regulation_l,'miscellaneous functions':other_l}
    return protFam

def protResults():
    catabolic_enzyme_l = [] ; isomerase_enzyme_l = [] ; binding_enzyme_l = [] ; ion_regulation_l = [] ; genetic_regulation_l = [] ; transcriptomic_regulation_l = []
    cytoskeleton_regulation_l = [] ; redox_regulation_l = [] ; immune_regulation_l = [] ; other_l = []
    protRes = {'catabolic enzyme':catabolic_enzyme_l,'isomerase enzyme':isomerase_enzyme_l,'binding enzyme':binding_enzyme_l,
    'ion regulation':ion_regulation_l,'genetic regulation':genetic_regulation_l,'translation regulation':transcriptomic_regulation_l,
    'cytoskeleton regulation':cytoskeleton_regulation_l,'redox regulation':redox_regulation_l,'immune regulation':immune_regulation_l,'miscellaneous functions':other_l}
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
        exprValue = 0
        if len(protExpr[i]) != 0:
            for j in range(len(protExpr[i])):
                exprValue += protExpr[i][j] 
        else:
            exprValue = 0
        protExpr[i] = exprValue 
    return protAnnot,protExpr 

def exploitResults(dfs,conditions,experiment):
    print(dfs)
    if len(conditions) == 0:
        print('No results are available for this experiment condition')
        return
    for i in range(len(dfs)):
        protAnnot,protExpr = sortResults(dfs[i],conditions[i])
        listProt = list(protExpr.keys()) 
        listExpr = list(protExpr.values()) 
        df = pd.DataFrame(list(zip(listProt, listExpr)),columns=['Proteins functions','Expression values'])
        
        df['colors'] = 'r'
        df.loc[df['Expression values']>=0,'colors'] = 'g'
        #print(df)
        
        fig = plt.figure()
        plt.bar(df['Proteins functions'], df['Expression values'], color=df['colors'],edgecolor='black')
        plt.axhline(y=0,color='black')
        plt.xticks(rotation=45,fontsize=10)
        fig.suptitle('Expression values of genes associated to proteins functions\n\n'+experiment+'_'+conditions[i], fontsize=20)
        plt.xlabel('Proteins functions', fontsize=16)
        plt.ylabel('Expression values\nadditive log2 fold change', fontsize=16)
        plt.subplots_adjust(bottom=0.3)
        #plt.show()
        

