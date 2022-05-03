"""
Final Expression Analyser

Marc Meynadier
"""
import os
import pandas as pd
import numpy as np

def retrieveResults():
    script_dir = os.path.dirname(__file__)
    full_path = os.path.join(script_dir, '../../data/net/8_functionnalAnnotation/functionnalGenesAnalysis/DESeq2_X_ontologizer')
    os.chdir(full_path)
    df = pd.read_csv('adult_preliminarySamples_pv_VS_gm_single_file_filtered.csv')
    return df

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

def sortResults():
    df = retrieveResults()
    protFam = protFamilies()
    protAnnot = protResults()
    protExpr = protResults()
    c=0
    for index, row in df.iterrows():   
        for j in protFam:
            for k in protFam[j]:
                if k in row['GO_annotation']:  
                    annot = row['pfam_annotation'].split(' ',1)[1]     
                    protAnnot[j].append(annot)   
                    protExpr[j].append(row['lfc_pv_VS_gm'])   
                    c+=1
                    break
    for i in protExpr:
        if len(protExpr[i]) != 0:
            averageExpr = np.mean(protExpr[i])
            protExpr[i] = averageExpr
    return protAnnot,protExpr
    
def exploitResults():
    protAnnot,protExpr = sortResults()
    print('\nProteins results :\n')
    print(protAnnot)
    print('\nExpression results :\n')
    print(protExpr,'\n')

exploitResults()