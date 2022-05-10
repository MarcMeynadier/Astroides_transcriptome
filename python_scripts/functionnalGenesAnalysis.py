"""
Functionnal genes analysis program

Three modules are available in this program :
1 : Retrieving DESeq2 results, sorting and concatenating them with the following information: Pfam annotation, protein sequences and Gene Ontology terms
2 : Ontologizer results retrieval, sorting and comparison across multiple files to detect common biological functions across multiple conditions
3 : Parsing and preparation of input files required to run the following functional enrichment programs: Ontologizer and GO_MWU


Marc Meynadier
"""

#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import sys
from DESeq2_analysis import *
from enrichment_analysis_parser import *
from ontologizer_analysis import *
from DESeq2_X_ontologizer import *
from getCandidates import *
from finalExpressionAnalyser import *


#------------------------------------------------------------------------------#
#                              Menu functions                                  #
#------------------------------------------------------------------------------#

def main_menu_display():
    print("\n")
    print("--------------------------------------------")
    print("|                                          |")
    print("|        Functionnal Genes Analysis        |")
    print("|                                          |")
    print("|                                          |")
    print("|       DESeq2_analyser : 1                |")
    print("|                                          |")
    print("|       Ontologizer_analyser: 2            |")
    print("|                                          |")
    print("|       DESeq2_X_ontologizer: 3            |")
    print("|                                          |")
    print("|       Final_Expression_analyser : 4      |") 
    print("|                                          |")
    print("|       Parsing & settings : 5             |")
    print("|                                          |")
    print("|       Exit : 6                           |")
    print("|                                          |")
    print("--------------------------------------------")
    print("\n")
    return

def menu_display_DESeq2():
    print("\n")
    print("----------------------------------------------")
    print("|                                            |")
    print("|        DESeq2 output files analysis        |")
    print("|                                            |")
    print("|                                            |")
    print("|       Single file annotation : 1           |")
    print("|                                            |")
    print("|       Genes unshared : 2                   |")
    print("|                                            |")
    print("|       Genes shared : 3                     |")
    print("|                                            |")
    print("|       Back to main menu : 4                |")
    print("|                                            |")
    print("----------------------------------------------")
    print("\n")
    return

def menu_display_ontologizer():
    print("\n")
    print("----------------------------------------------")
    print("|                                            |")
    print("|     Ontologizer output files analysis      |")
    print("|                                            |")
    print("|                                            |")
    print("|       Single file enrichment : 1           |")
    print("|                                            |")
    print("|       GO terms unshared : 2                |")
    print("|                                            |")
    print("|       GO terms shared : 3                  |")
    print("|                                            |")
    print("|       Back to main menu : 4                |")
    print("|                                            |")
    print("----------------------------------------------")
    print("\n")
    return

def menu_display_parsing_settings():
    print("\n")
    print("----------------------------------------------")
    print("|                                            |")
    print("|   Parsing input files & program settings   |")
    print("|                                            |")
    print("|                                            |")
    print("|       Ontologizer input files parsing : 1  |")
    print("|                                            |")
    print("|       Generate candidate genes file : 2    |")
    print("|                                            |")
    print("|       p-value threshold : 3                |")
    print("|                                            |")
    print("|       Candidate genes threshold : 4        |")
    print("|                                            |")
    print("|       Back to main menu : 5                |")
    print("|                                            |")
    print("----------------------------------------------")
    print("\n")
    return

def menu_display_enrichment_parsing():
    print("\n")
    print("----------------------------------------------")
    print("|                                            |")
    print("|     Enrichment analysis files parsing      |")
    print("|                                            |")
    print("|                                            |")
    print("|       Ontologizer : 1                      |")
    print("|                                            |")
    print("|       GO_MWU : 2                           |")
    print("|                                            |")
    print("|       Back to main menu : 3                |")
    print("|                                            |")
    print("----------------------------------------------")
    print("\n")
    return

default_threshold = 0.05
defaultFlagCandidate = 'N'

def main_menu(threshold,flagCandidate):
    while True:
        scriptDir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(scriptDir) 
        main_menu_display()
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 6\n")
                continue
            break 
        if answer==1:
            menu_DESeq2(threshold,flagCandidate)
        elif answer==2:
            menu_ontologizer_output(threshold,flagCandidate)
        elif answer==3:
            matchingFiles() 
        elif answer==4:
            typeOrg, experiment, org= experimentChoice() 
            filenames = getFilenamesFinal(experiment,threshold,flagCandidate)
            if len(filenames)==0:
                print('\nNo results are available for this experiment condition')
                main_menu(threshold,flagCandidate)
            dfs,conditions,experiment = filenamesToDfFinal(filenames,experiment,flagCandidate)
            exploitResults(dfs,conditions,experiment)
        elif answer==5:
            parsing_settings(threshold,flagCandidate)
        elif answer==6:
            sys.exit(0)

def menu_DESeq2(threshold,flagCandidate):
    typeOrg, experiment, org = experimentChoice()
    filenames = getFilenames(typeOrg,experiment)
    while True:
        menu_display_DESeq2()
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 4\n")
                continue
            break
        if answer==1:
            singleFile(filenames,experiment,org,threshold,flagCandidate)
        elif answer==2:
            genesUnshared(filenames,experiment,org,threshold,flagCandidate)
        elif answer==3:
            genesShared(filenames,experiment,org,threshold,flagCandidate)
        elif answer==4:
            main_menu(threshold,flagCandidate)

def menu_ontologizer_output(threshold,flagCandidate):
    while True:
        menu_display_ontologizer() 
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 4\n")
                continue
            break
        scriptDir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(scriptDir) 
        if answer==1:
            singleFileOntologizer(threshold,flagCandidate)
        elif answer==2:
            genesUnsharedOntologizer(threshold,flagCandidate)
        elif answer==3:
            genesSharedOntologizer(threshold,flagCandidate)
        elif answer==4:
            main_menu(threshold,flagCandidate)


def parsing_settings(threshold,flagCandidate):
        while True: 
            menu_display_parsing_settings()
            while True:
                try:
                    answer = int(input())
                except ValueError:
                    print("\nYou must indicate an integer value ranging from 1 to 5\n")
                    continue
                break
            if answer==1:
                menu_enrichment_parsing(threshold,flagCandidate)
            elif answer==2:
                getCandidateGenes()
            elif answer==3:
                threshold = setThreshold()
            elif answer==4:
                flagCandidate = filterByCandidate()
            elif answer==5:
                main_menu(threshold,flagCandidate) 

def menu_enrichment_parsing(threshold,flagCandidate):
    while True:
        menu_display_enrichment_parsing()
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 3\n")
                continue
            break
        if answer==1:
            getOntologyFileOntologizer()
            getAssociationFile()
            getPopulationFile()
            getStudysetFileOntologizer(threshold)
        elif answer==2:
            getOntologyFileGOMWU()
            getStudysetFileGOMWU(threshold)
        elif answer==3:
            main_menu(threshold,flagCandidate)


#------------------------------------------------------------------------------#
#                                    Main                                      #
#------------------------------------------------------------------------------#


main_menu(default_threshold,defaultFlagCandidate)
