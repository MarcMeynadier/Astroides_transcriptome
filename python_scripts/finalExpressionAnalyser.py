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
import statistics
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sb
from matplotlib.backends.backend_pdf import PdfPages
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
    path = os.path.join(script_dir, '../../data/net/8_functionalAnnotation/functionalGenesAnalysis/DESeq2_X_ontologizer')
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
    comparedFiles
        list, contains the Pandas dataframes for each selected files.
    conditions,
        list, contains names of the conditions of experimentation of each selected files.
    experiment
        str, contains the type of experiment retrieved with experimentChoice(). 
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
    """
    Description
    -----------
    Creates a dictionary whose keys are biological functions and the values of the keywords associated with these same functions. 
    Each keyword has been carefully selected from the Gene Ontology terms, and its membership to a biological function is documented. 

    Parameters
    ----------
    None

    Returns
    -------
    protFam
        dictionnary, each key corresponds to a biological function, and each value to a keyword of the Gene ontology.
    """

    catabolic_enzyme_l=['hydrolase','hydrolysis','protease','peptidase','endopeptidase','metalloendopeptidase','proteolysis','phosphatase','phospholipase'
    ,'catalytic activity','metallopeptidase','metallocarboxypeptidase','lipase','catabolic process','ATPase','GTPase','aminidase']
    isomerase_enzyme_l=['lyase','isomerase','transferase','mutase','oxydoreductase','dehydrogenase','kinase','peroxydase','deiodinase']
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
    """
    Description
    -----------
    Create a dictionary identical to protFamilies() but whose keys do not contain any values.  

    Parameters
    ----------
    None

    Returns
    -------
    protRes
        dictionnary, each key corresponds to a biological function, and no values are added to the keys.
    """

    catabolic_enzyme_l = [] ; isomerase_enzyme_l = [] ; binding_enzyme_l = [] ; ion_regulation_l = [] ; genetic_regulation_l = [] ; transcriptomic_regulation_l = []
    cytoskeleton_regulation_l = [] ; redox_regulation_l = [] ; immune_regulation_l = [] ; other_l = []
    protRes = {'catabolic enzyme':catabolic_enzyme_l,'isomerase enzyme':isomerase_enzyme_l,'binding enzyme':binding_enzyme_l,
    'ion regulation':ion_regulation_l,'genetic regulation':genetic_regulation_l,'translation regulation':transcriptomic_regulation_l,
    'cytoskeleton regulation':cytoskeleton_regulation_l,'redox regulation':redox_regulation_l,'immune regulation':immune_regulation_l,'miscellaneous functions':other_l}
    return protRes

def sortResultsBarplots(df,condition):
    """
    Description
    -----------
    From the filled keyword dictionary provided by protaFamilies(), the empty keyword dictionary provided by protResults() and the list of dataframes 
    of each file provided by filenamesToDfFinal(), sort the Gene Ontology terms of each dataframe in order to match a gene to its biological function(s). 
    Two empty dictionaries are created with protResults(), one allowing to store the proteins corresponding to the functions, the other allowing to store 
    the expression values of the genes. Once all the dataframes have been scanned, the two dictionaries are returned. 

    Parameters
    ----------
    df
        Pandas dataframe, contains the information of an output file of DESeq2_X_ontologizer.
    condition
        str, contains name of the conditions of experimentation of the selected files.

    Returns
    -------
    protAnnot
        dictionnary, contains the names of the proteins associated with each biological function for a DESeq2_X_ontologizer output file.
    protExpr
        dictionnary, contains the genes expressions values associated with each biological function for a DESeq2_X_ontologizer output file. 
    """

    protFam = protFamilies()
    protAnnot = protResults()
    protExpr = protResults()
    exprSD = []
    for index, row in df.iterrows():   
        for j in protFam:
            for k in protFam[j]:
                if k in row['GO_annotation']:  
                    annot = row['pfam_annotation'].split(' ',1)[1]     
                    protAnnot[j].append(annot)   
                    protExpr[j].append(row['lfc_'+condition])   
                    break
    for i in protExpr:
        exprValue = 0
        exprList = []
        if len(protExpr[i]) != 0:
            for j in range(len(protExpr[i])):
                exprList.append(protExpr[i][j])
                exprValue += protExpr[i][j]  
        else:
            exprValue = 0
            exprList.append(exprValue)
        exprSD.append(np.std(exprList))
        protExpr[i] = exprValue 
    return protAnnot,protExpr,exprSD

def sortResultsBoxplots(df,condition):
    protFam = protFamilies()
    protAnnot = protResults()
    protExpr = protResults()
    df = df.rename(columns={df.columns[3]: 'lfc'})
    #df = df[df.lfc < 10]
    for index, row in df.iterrows():   
        for j in protFam:
            for k in protFam[j]:
                if k in row['GO_annotation']:  
                    annot = row['pfam_annotation'].split(' ',1)[1]     
                    protAnnot[j].append(annot)   
                    protExpr[j].append(row['lfc'])   
                    break
    return protAnnot,protExpr

def exploitResultsBarplots(dfs,conditions,experiment):
    """
    Description
    -----------
    Using the Matplotlib library, create a bar chart representing the different biological functions of each file according to the expression values 
    of each biological function. It is possible to visualize on the bar chart up to 4 different files, a type of pattern represented in legend being 
    displayed for each file. The bar chart thus created is saved then displayed. 

    Parameters
    ----------
    dfs
        list, contains the Pandas dataframes for each selected files.
    conditions,
        list, contains names of the conditions of experimentation of each selected files.
    experiment
        str, contains the type of experiment retrieved with experimentChoice().
        
    Returns
    -------
    None
    """

    if len(conditions) == 0:
        print('No results are available for this experiment condition')
        return
    outputPath = '../../../../../output/functionalGenesAnnotation/'
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
        protAnnot,protExpr,exprSD = sortResultsBarplots(dfs[i],conditions[i])
        protDf = pd.DataFrame.from_dict(protAnnot,orient='index') 
        protDf2 = pd.DataFrame.from_dict(protExpr,orient='index')
        protDf2.columns=['lfc'] ; protDf2['lfc'].astype(str)
        protDf = pd.concat([protDf,protDf2],axis=1) 
        protDf = protDf.transpose() 
        protDf.to_csv(outputPath+'proteinsTableAdditiveLFC_'+experiment+'_'+conditions[i]+'.csv')
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
        df.insert(3,"Standard deviation of expression",exprSD,True)
        bars = plt.bar(x_axis+tick, df['Expression values '+conditions[i]], color=df['colors '+conditions[i]],edgecolor='black',width=width,hatch=patterns[i])
        plt.errorbar(x_axis+tick,df['Expression values '+conditions[i]], xerr=None,yerr=df['Standard deviation of expression'], ecolor='black', color='black',ls='none',elinewidth=2)
        bars_list.append(bars)
    axes = plt.gca() ; ymax = axes.get_ylim()  
    if ymax[1] < 0:
       y_neg_coeff = 2 
    elif ymax[1] < 10 and ymax[0] < -10:
       y_coeff = 2 ; y_neg_coeff = 4
    elif ymax[1] <= 5:
        y_coeff = 0.15 ; y_neg_coeff = 0.8
    elif ymax[1] > 5 and ymax[1] <= 10:
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
    elif ymax[1] > 110 and ymax[1] <= 135:
        y_coeff = 150 ; y_neg_coeff = 400   
    textBarValue = ymax[1] 
    if textBarValue == 0:
        textBarValue = 8
    upLegendSwitch = False
    downLegendSwitch = False
    for i in range(len(bars_list)):
        for j in range(len(bars_list[i])):
                yval = bars_list[i][j].get_height() 
                if yval >= 0:
                    upLegendSwitch = True
                    plt.text(bars_list[i][j].get_x()+adjust, yval+(y_coeff*(1/textBarValue)), 'N:'+str(numberProtList[i][j]),ha='center',fontsize=fontsize)
                else:
                    downLegendSwitch = True
                    plt.text(bars_list[i][j].get_x()+adjust,yval-(y_neg_coeff*(1/textBarValue)), 'N:'+str(numberProtList[i][j]),ha='center',fontsize=fontsize)
    if len(dfs) == 1:
        plt.xticks(x_axis+0.3,prot_functions,rotation=45,fontsize=12) 
    elif len(dfs) == 2:
        plt.xticks(x_axis+0.6-(0.3/len(dfs)),prot_functions,rotation=45,fontsize=12)     
    elif len(dfs) == 3:
        plt.xticks(x_axis+(0.3*len(dfs)-0.5),prot_functions,rotation=45,fontsize=12)
    elif len(dfs) == 4:     
        plt.xticks(x_axis+(0.2*len(dfs)-0.3),prot_functions,rotation=45,fontsize=12)  
    plt.axhline(y=0,color='black')
    plt.xlabel('Protein functions', fontsize=18)
    plt.ylabel('Additive log₂ fold change', fontsize=18)
    legend = []
    if upLegendSwitch == True:
        up = mpatches.Patch(facecolor='#006CD1', label='Upregulated',edgecolor='black') ; legend.append(up) # juvenile
    if downLegendSwitch == True:
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
        elif i==3:
            leg = mpatches.Patch(hatch='***',label=conditions[i],facecolor='white',edgecolor='black')
            legend.append(leg)     
    plt.legend(handles=legend, loc=4,fontsize=14,frameon=False)
    #plt.title('Expression values of genes associated to proteins functions\n\n'+experiment,fontsize=18) # juvenile
    plt.subplots_adjust(bottom=0.3) 
    conditionsStr = '_X_'.join(conditions) 
    fig.savefig(outputPath+'FGA_barplot_'+experiment+'_'+conditionsStr+'.png')   
    plt.show()

def prepareBoxplots(dfs,conditions):
    catabolic_enzyme_df = pd.DataFrame() ; isomerase_enzyme_df = pd.DataFrame() ; binding_enzyme_df = pd.DataFrame() 
    ion_regulation_df = pd.DataFrame() ; genetic_regulation_df = pd.DataFrame() ; transcriptomic_regulation_df = pd.DataFrame() 
    cytoskeleton_regulation_df = pd.DataFrame() ; redox_regulation_df = pd.DataFrame() ; immune_regulation_df = pd.DataFrame()  
    other_df = pd.DataFrame() 
    for i in range(len(dfs)):
        protAnnot,protExpr = sortResultsBoxplots(dfs[i],conditions[i])
        #protDf = pd.DataFrame.from_dict(protAnnot,orient='index') 
        exprDf = pd.DataFrame.from_dict(protExpr,orient='index') 
        catabolic_enzyme_df[conditions[i]] = exprDf.loc['catabolic enzyme']
        isomerase_enzyme_df[conditions[i]] = exprDf.loc['isomerase enzyme']
        binding_enzyme_df[conditions[i]] = exprDf.loc['binding enzyme']
        ion_regulation_df[conditions[i]] = exprDf.loc['ion regulation']
        genetic_regulation_df[conditions[i]] = exprDf.loc['genetic regulation']
        transcriptomic_regulation_df[conditions[i]] = exprDf.loc['translation regulation']
        cytoskeleton_regulation_df[conditions[i]] = exprDf.loc['cytoskeleton regulation']
        redox_regulation_df[conditions[i]] = exprDf.loc['redox regulation']
        immune_regulation_df[conditions[i]] = exprDf.loc['immune regulation']
        other_df[conditions[i]] = exprDf.loc['miscellaneous functions']
    catabolic_enzyme_df = catabolic_enzyme_df.assign(Function="Catabolic enzyme")
    isomerase_enzyme_df = isomerase_enzyme_df.assign(Function="Isomerase enzyme")
    binding_enzyme_df = binding_enzyme_df.assign(Function="Binding enzyme")
    ion_regulation_df = ion_regulation_df.assign(Function="Ionic regulation")
    genetic_regulation_df = genetic_regulation_df.assign(Function="Genetic regulation")
    transcriptomic_regulation_df = transcriptomic_regulation_df.assign(Function="Transcriptomic regulation")
    cytoskeleton_regulation_df = cytoskeleton_regulation_df.assign(Function="Cytoskeleton regulation")
    redox_regulation_df = redox_regulation_df.assign(Function="Redox regulation")
    other_df = other_df.assign(Function="Miscellaneous functions")
    concatDf = pd.concat([catabolic_enzyme_df,isomerase_enzyme_df,binding_enzyme_df,ion_regulation_df,
    genetic_regulation_df,transcriptomic_regulation_df,cytoskeleton_regulation_df,redox_regulation_df,other_df]) 
    meltDf = pd.melt(concatDf,id_vars=['Function'],var_name=['Conditions'])
    return meltDf

def exploitResultsBoxplots(meltDf,conditions,experiment):
    if len(conditions) == 0:
        print('No results are available for this experiment condition')
        return 
    plt.figure(figsize=(12,10))
    ax = sb.boxplot(x="Function", y="value", hue="Conditions",data=meltDf)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 6.5)
    plt.setp(ax.get_xticklabels(), rotation=30)
    plt.xlabel("Protein functions") ; plt.ylabel("LFC")
    conditionsStr = '_X_'.join(conditions) 
    outputPath = '../../../../../output/functionalGenesAnnotation/' 
    plt.savefig(outputPath+'FGA_boxplot_'+experiment+'_'+conditionsStr+'.png')
    plt.show() 
    plt.close()
    
    
    