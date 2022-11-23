"""
Functional genes analysis program

Five modules are available in this program :
1 : Retrieving DESeq2 results, sorting and concatenating them with the following information: Pfam annotation, protein sequences and Gene Ontology terms
2 : Ontologizer results retrieval, sorting and comparison across multiple files to detect common biological functions across multiple conditions
3 : Parsing and preparation of input files required to run the following functional enrichment programs: Ontologizer and GO_MWU

Marc Meynadier
"""

#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import sys
from DESeq2Analysis import * 
from enrichmentAnalysisParser import *
from ontologizerAnalysis import *
from DESeq2_X_ontologizer import *
from thresholdSettings import *
from getCandidates import *
from finalExpressionAnalyser import *


#------------------------------------------------------------------------------#
#                              Menu functions                                  #
#------------------------------------------------------------------------------#


def main_menu_display():
    """
    Description
    -----------
    Shows the user the main menu of the program. 

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    print("\n")
    print("--------------------------------------------")
    print("|                                          |")
    print("|        Functional Genes Analysis         |")
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
    """
    Description
    -----------
    Shows the user the menu of the DESeq2Analyser module. 

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

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
    """
    Description
    -----------
    Shows the user the menu of the ontologizerAnalyser module. 

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

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
    """
    Description
    -----------
    Shows the user the menu of settings module. 

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

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
    """
    Description
    -----------
    Shows the user the menu of parsing settings. 

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

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

def main_menu(threshold,flagCandidate):
    """
    Description
    -----------
    Main menu of the program to call the other modules. 

    Parameters
    ----------
    threshold
        int, contains the value of p-value threshold defined in the settings. A default value of 0.05 is set.
    flagCandidate
        str, contains the switch activating the filtering by candidate genes defined in the settings. A default value of N (no) is set.

    Returns
    -------
    None
    """

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
            filenames = getFilenamesFinal(experiment,flagCandidate)
            if len(filenames)==0:
                print('\nNo results are available for this experiment condition')
                main_menu(threshold,flagCandidate)
            dfs,conditions,experiment = filenamesToDfFinal(filenames,experiment,flagCandidate)
            prepareBoxplots(dfs,conditions)
            #exploitResultsBoxplots(dfs,conditions,experiment)
        elif answer==5:
            parsing_settings(threshold,flagCandidate)
        elif answer==6:
            sys.exit(0)

def menu_DESeq2(threshold,flagCandidate):
    """
    Description
    -----------
    Main menu of the DESeq2Analyser module. 

    Parameters
    ----------
    threshold
        int, contains the value of p-value threshold defined in the settings. 
    flagCandidate
        str, contains the switch activating the filtering by candidate genes defined in the settings. 

    Returns
    -------
    None
    """

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
    """
    Description
    -----------
    Main menu of the ontologizerAnalyser module. 

    Parameters
    ----------
    threshold
        int, contains the value of p-value threshold defined in the settings. 
    flagCandidate
        str, contains the switch activating the filtering by candidate genes defined in the settings. 

    Returns
    -------
    None
    """

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
    """
    Description
    -----------
    Main menu of settings module. 

    Parameters
    ----------
    threshold
        int, contains the value of p-value threshold defined in the settings. 
    flagCandidate
        str, contains the switch activating the filtering by candidate genes defined in the settings. 

    Returns
    -------
    None
    """

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
    """
    Description
    -----------
    Main menu of parsing settings module. 

    Parameters
    ----------
    threshold
        int, contains the value of p-value threshold defined in the settings. 
    flagCandidate
        str, contains the switch activating the filtering by candidate genes defined in the settings. 

    Returns
    -------
    None
    """

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


default_threshold = 0.05
defaultFlagCandidate = 'N'
main_menu(default_threshold,defaultFlagCandidate)
