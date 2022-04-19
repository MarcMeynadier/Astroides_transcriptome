"""
Comparison of overlaping genes betwen Ontologizer output files
Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import os
import sys
from regex import F
import wget
import pandas as pd
import csv
from urllib.request import urlopen


def getFilesNames():
    threshold_pvalue=0.05
    os.chdir('../data/net/7_functionnalAnnotation/ontologizer/outputResults/codingDEG') # Changing working directory to DESeq2 results
    path=os.getcwd()
    for file in os.listdir(path):
        if file.endswith(".txt"):
            curedFile = []
            fileName = file.replace('.csv','')
            csvFile = pd.read_csv(file,on_bad_lines='skip',sep='\t')
            print(csvFile)
            

getFilesNames()