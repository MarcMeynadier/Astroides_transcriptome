"""
Matching between differentially expressed genes retrieved from DESeq2
Marc Meynadier
"""
import os
import glob
import pandas as pd
from torch import threshold

os.chdir('../data/net/6_deseq2/adult') # Changing working directory to DESeq2 results
path=os.getcwd()
filenamesPreSamp = glob.glob(path + "/*preliminarySamples*.csv")
filenamesTrueTransp = glob.glob(path + "/*trueTransplant*.csv")

threshold=0.05

dfsPreSamp = [pd.read_csv(filename) for filename in filenamesPreSamp]
for i in range(len(dfsPreSamp)):
    dfsPreSamp[i]=dfsPreSamp[i].drop(dfsPreSamp[i][dfsPreSamp[i].padj>threshold].index)
    dfsPreSamp[i]=dfsPreSamp[i][dfsPreSamp[i]['padj'].notna()]
concatPreSamp = pd.concat(dfsPreSamp, ignore_index=True)
print(concatPreSamp)

dfsTrueTransp = [pd.read_csv(filename) for filename in filenamesTrueTransp]
for i in range(len(dfsTrueTransp)):
    dfsTrueTransp[i]=dfsTrueTransp[i].drop(dfsTrueTransp[i][dfsTrueTransp[i].padj>threshold].index)
    dfsTrueTransp[i]=dfsTrueTransp[i][dfsTrueTransp[i]['padj'].notna()]
concatTrueTransp = pd.concat(dfsTrueTransp, ignore_index=True)
print(concatTrueTransp)