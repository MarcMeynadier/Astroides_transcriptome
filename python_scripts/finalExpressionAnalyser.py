"""
finalExpressionAnalyser : Classifies DESeq2_X_ontologizer genes into different associated 
biological functions and produces barplots of these biological functions based on 
the gene expression of the affected genes. 

Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


#------------------------------------------------------------------------------#
#                              Files management                                #
#------------------------------------------------------------------------------#


def getFilenamesFinal(experiment,flagCandidate):
    """
    Description
    -----------
    Retrieves the output file names from DESeq2_X_ontologizer based on associated experiment.

    Parameters
    ----------
    experiment
        str, contains the type of experiment retrieved with experimentChoice().      
    flagCandidate
        str, contains the switch activating the filtering by candidate genes defined in the settings.

    Returns
    -------
    filenames
        list, contains the file names of the DESeq2_X_ontologizer results.
    """   

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
    """
    Description
    -----------
    Retrieves CSV files from DESeq2_X_ontologizer by their names, transforms them into dataframe before entering them into a list. 
    The genes are then filtered according to the p-value and candidate gene thresholds, and the dataframe list is returned.

    Parameters
    ----------
    filenames
        list, contains the name of files parsed by getFilenamesFinal().
    experiment
        str, contains the type of experiment retrieved with experimentChoice().      
    flagCandidate
        str, contains the switch activating the filtering by candidate genes defined in the settings.

    Returns
    -------
    filenames
        list, contains the file names of the DESeq2_X_ontologizer results.
    """   

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
    print("\nHow many files do you want to analyse ?\n")  
    nbrFiles = 3
    while nbrFiles!=0:
        try:
            nbrFiles = int(input())
            if nbrFiles == 0:
                print("\nYou must chose at least one file\n")
            else:
                break
        except ValueError:
            print("\nYou must indicate an integer number corresponding to the number of files you want to analyse\n")
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


#------------------------------------------------------------------------------#
#                            Analysis computation                              #
#------------------------------------------------------------------------------#


def protFamilies():
    catabolic_enzyme_l=['hydrolase','hydrolysis','protease','peptidase','endopeptidase','metalloendopeptidase','proteolysis','phosphatase','phospholipase'
    ,'catalytic activity','metallopeptidase','metallocarboxypeptidase','lipase','catabolic process','ATPase','GTPase','aminidase']
    isomerase_enzyme_l=['lyase','isomerase','transferase','mutase','oxydoreductase','dehydrogenase','kinase,peroxydase,deiodinase']
    binding_enzyme_l=['carbohydrate binding','chitin binding','ion binding','acid binding','protein binding','biotin binding']
    ion_regulation_l=['proton','ion transport','calcium ion','zinc ion','magnesium ion','channel','transmembrane transport','intermembrane','translocase',
    'iron ion','copper ion','cation transport','anion transport','sulfate transport','potassium channel','ammonium transmembrane','ionotropic glutamate receptor',
    'ion transmembrane']
    genetic_regulation_l=['DNA','histone','nucleic acid','telomere','helicase']
    transcriptomic_regulation_l=['translation','ribosome','RNA','elongation']
    cytoskeleton_regulation_l = ['microtubule','cytoskeleton','mitotic spindle','cytokinesis','extracellular matrix','actin','cell adhesion','extracellular matrix','dynein'
    ,'dystrophin','fibrillin','sushi','galaxin']
    redox_regulation_l = ['oxydoreductase','oxygenase','antioxydant','redox homeostasis','NAD','NADH','dehydrogenase','FAD','FADH','peroxydase','deiodinase']
    immune_regulation_l = ['tumor','necrosis','apoptosis','immune']
    other_l = ['signal transduction','vesicle','lipid transport','organism development','nuclear pore','circadian rhythm']
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
    if len(conditions) == 0:
        print('No results are available for this experiment condition')
        return
    tick = 0 ; fontsize = 0 ; width = 0 ; adjust = 0
    x_axis = np.array
    prot_functions = pd.Series
    patterns = ["//", "..", "--", "**" ]
    fig = plt.figure(figsize=(15, 10))
    bars_list = []
    numberProtList = []
    for i in range(len(dfs)):
        if len(dfs) < 3:
            tick += 0.3
            width = 0.3
            fontsize = 8
            adjust = 0.15
        else:
            width = 0.2
            tick += 0.2
            fontsize = 6
            adjust = 0.1
        protAnnot,protExpr = sortResults(dfs[i],conditions[i])
        print(protAnnot) ; print(len(protAnnot))
        numberProt = []
        for j in protAnnot.values():
            numberProt.append(len(j))
        numberProtList.append(numberProt)
        listProt = list(protExpr.keys()) 
        listExpr = list(protExpr.values()) 
        df = pd.DataFrame(list(zip(listProt, listExpr)),columns=['Proteins functions','Expression values '+conditions[i]])  
        df['colors '+conditions[i]] = '#E66100'
        df.loc[df['Expression values '+conditions[i]]>=0,'colors '+conditions[i]] = '#006CD1'       
        x_axis = np.arange(len(df['Proteins functions']))
        prot_functions = df['Proteins functions']
        bars = plt.bar(x_axis+tick, df['Expression values '+conditions[i]], color=df['colors '+conditions[i]],edgecolor='black',width=width,hatch=patterns[i])
        bars_list.append(bars)
    axes = plt.gca() ; ymax = axes.get_ylim()  ; print(ymax)
    if ymax[1] <= 5:
        y_coeff = 0.15 ; y_neg_coeff = 0.8
    if ymax[1] > 5 and ymax[1] <= 10:
        y_coeff = 0.5 ; y_neg_coeff = 1.5
    elif ymax[1] > 10 and ymax[1] <= 30:
        y_coeff = 2 ; y_neg_coeff = 6
    elif ymax[1] > 30 and ymax[1] <= 50: 
        y_coeff = 20 ; y_neg_coeff = 80
    elif ymax[1] > 50 and ymax[1] <= 70:
        y_coeff = 60 ; y_neg_coeff = 120
    elif ymax[1] > 70 and ymax[1] <= 90:
        y_coeff = 90 ; y_neg_coeff = 270
    elif ymax[1] > 90 and ymax[1] <= 110:
        y_coeff = 120 ; y_neg_coeff = 350
    elif ymax[1] > 110 and ymax[1] <= 130:
        y_coeff = 150 ; y_neg_coeff = 400   
    for i in range(len(bars_list)):
        for j in range(len(bars_list[i])):
                yval = bars_list[i][j].get_height()
                if yval >= 0:
                    plt.text(bars_list[i][j].get_x()+adjust, yval+(y_coeff*(1/ymax[1])), 'N:'+str(numberProtList[i][j]),ha='center',fontsize=fontsize)
                else:
                    plt.text(bars_list[i][j].get_x()+adjust,yval-(y_neg_coeff*(1/ymax[1])), 'N:'+str(numberProtList[i][j]),ha='center',fontsize=fontsize)
    if len(dfs) == 1:
        plt.xticks(x_axis+0.3,prot_functions,rotation=45,fontsize=10) 
    elif len(dfs) == 2:
        plt.xticks(x_axis+0.6-(0.3/len(dfs)),prot_functions,rotation=45,fontsize=10)     
    elif len(dfs) == 3:
        plt.xticks(x_axis+(0.3*len(dfs)-0.3),prot_functions,rotation=45,fontsize=10)
    elif len(dfs) == 4:     
        plt.xticks(x_axis+(0.2*len(dfs)-0.3),prot_functions,rotation=45,fontsize=10)  
    plt.axhline(y=0,color='black')
    plt.xlabel('Proteins functions', fontsize=16)
    plt.ylabel('Expression values\n(additive log2 fold change)', fontsize=16)
    legend = []
    up = mpatches.Patch(facecolor='#006CD1', label='Upregulated',edgecolor='black') ; legend.append(up)
    down = mpatches.Patch(facecolor='#E66100', label='Downregulated',edgecolor='black') ; legend.append(down)
    for i in range(len(conditions)):
        if i==0:
            leg = mpatches.Patch(hatch=r'////',label=conditions[i],facecolor='white',edgecolor='black')
            legend.append(leg)
        elif i==1:
            leg = mpatches.Patch(hatch='....',label=conditions[i],facecolor='white',edgecolor='black')
            legend.append(leg)
        elif i==2:
            leg = mpatches.Patch(hatch='----',label=conditions[i],facecolor='white',edgecolor='black')
            legend.append(leg)
        elif i==2:
            leg = mpatches.Patch(hatch='****',label=conditions[i],facecolor='white',edgecolor='black')
            legend.append(leg)     
    plt.legend(handles=legend, loc=0,fontsize=7.5,frameon=False)
    plt.title('Expression values of genes associated to proteins functions\n\n'+experiment,fontsize=18)
    plt.subplots_adjust(bottom=0.3)
    outputPath = '../../../../../output/functionalGenesAnnotation/'
    conditionsStr = '_X_'.join(conditions) 
    fig.savefig(outputPath+'FGA_barplot_'+experiment+'_'+conditionsStr+'.png')   
    #plt.show()
        
