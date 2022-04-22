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


#------------------------------------------------------------------------------#
#                              Menu functions                                  #
#------------------------------------------------------------------------------#

def main_menu_display():
    print("\n")
    print("--------------------------------------------------------")
    print("|                                                      |")
    print("|             Functionnal Genes Analysis               |")
    print("|                                                      |")
    print("|                                                      |")
    print("|      DESeq2 output files analysis : 1                |")
    print("|                                                      |")
    print("|      Ontologizer output files analysis : 2           |")
    print("|                                                      |")
    print("|      DESeq2 output filtering after enrichment : 3    |")
    print("|                                                      |")
    print("|      Enrichment analysis input files parsing : 4     |") 
    print("|                                                      |")
    print("|      Exit : 5                                        |")
    print("|                                                      |")
    print("--------------------------------------------------------")
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
    print("|       Back to previous selection : 4       |")
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
    print("|       Back to previous selection : 4       |")
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
    print("|       Ontologizer  : 1                     |")
    print("|                                            |")
    print("|       GO_MWU : 2                           |")
    print("|                                            |")
    print("|       Back to previous selection : 3       |")
    print("|                                            |")
    print("----------------------------------------------")
    print("\n")
    return

def main_menu():
    while True:
        scriptDir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(scriptDir) 
        main_menu_display()
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 5\n")
                continue
            break
        if answer==1:
            menu_DESeq2()
        elif answer==2:
            menu_ontologizer_output()
        elif answer==3:
            matchingFiles()
        elif answer==4:
            menu_enrichment_parsing()
        elif answer==5:
            sys.exit(0)

def menu_DESeq2():
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
            singleFile(filenames,experiment,org)
        elif answer==2:
            genesUnshared(filenames,experiment,org)
        elif answer==3:
            genesShared(filenames,experiment,org)
        elif answer==4:
            main_menu()

def menu_ontologizer_output():
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
            singleFileOntologizer()
        elif answer==2:
            genesUnsharedOntologizer()
        elif answer==3:
            genesSharedOntologizer()
        elif answer==4:
            main_menu()

def menu_enrichment_parsing():
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
            getStudysetFileOntologizer()
        elif answer==2:
            getOntologyFileGOMWU()
            getStudysetFileGOMWU()
        elif answer==3:
            main_menu()


#------------------------------------------------------------------------------#
#                                    Main                                      #
#------------------------------------------------------------------------------#


main_menu()
