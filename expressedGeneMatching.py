"""
Matching between differentially expressed genes retrieved from DESeq2
Marc Meynadier
"""
import os
import glob
import pandas as pd
import numpy as np

threshold_pvalue=0.05

os.chdir('../data/net/6_deseq2/adult') # Changing working directory to DESeq2 results
path=os.getcwd()
filenamesPreSamp = glob.glob(path + "/*preliminarySamples*.csv")
filenamesTrueTransp = glob.glob(path + "/*trueTransplant*.csv")

threshold=0.05

dfsPreSamp = [pd.read_csv(filename) for filename in filenamesPreSamp]
for i in range(len(dfsPreSamp)):
    dfsPreSamp[i]=dfsPreSamp[i].drop(dfsPreSamp[i][dfsPreSamp[i].padj>threshold_pvalue].index)
    dfsPreSamp[i]=dfsPreSamp[i][dfsPreSamp[i]['padj'].notna()]
concatPreSamp = pd.concat(dfsPreSamp, ignore_index=True)

dfsTrueTransp = [pd.read_csv(filename) for filename in filenamesTrueTransp]
for i in range(len(dfsTrueTransp)):
    dfsTrueTransp[i]=dfsTrueTransp[i].drop(dfsTrueTransp[i][dfsTrueTransp[i].padj>threshold_pvalue].index)
    dfsTrueTransp[i]=dfsTrueTransp[i][dfsTrueTransp[i]['padj'].notna()]
concatTrueTransp = pd.concat(dfsTrueTransp, ignore_index=True)
#print(concatPreSamp)
#print(concatTrueTransp)

#for index, row in concatPreSamp.iterrows():
#  print(row[0])


def comparingDataframe(df1,df2):
    sharedGenes=[]
    unsharedGenes=[]
    for i, row in df1.iterrows():
        for j, row2 in df2.iterrows():
            if row[0]==row2[0]:
                sharedGenes.append(row[0])
    sharedGenes=np.unique(sharedGenes)
    print("list of shared genes :",sharedGenes)
    print("Number of shared genes :",len(sharedGenes))
    return sharedGenes

            
comparingDataframe(concatPreSamp,concatTrueTransp)