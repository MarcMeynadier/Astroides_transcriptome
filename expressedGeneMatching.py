"""
Matching between differentially expressed genes retrieved from DESeq2
Marc Meynadier
"""
import os
import glob
import pandas as pd
import numpy as np
from functools import reduce

threshold_pvalue=0.05

os.chdir('../data/net/6_deseq2/adult') # Changing working directory to DESeq2 results
path=os.getcwd()
filenamesTrueTransp = glob.glob(path + "/*trueTransplant*.csv")

for filename in filenamesTrueTransp:
    print(filename)

dfsTrueTransp = [pd.read_csv(filename) for filename in filenamesTrueTransp]
for i in range(len(dfsTrueTransp)):
    dfsTrueTransp[i].rename(columns={ dfsTrueTransp[i].columns[0]: "gene" }, inplace = True)
    dfsTrueTransp[i]=dfsTrueTransp[i].drop(dfsTrueTransp[i][dfsTrueTransp[i].padj>threshold_pvalue].index)
    dfsTrueTransp[i]=dfsTrueTransp[i][dfsTrueTransp[i]['padj'].notna()]
concatTrueTransp = pd.concat(dfsTrueTransp, ignore_index=True)
#print(concatPreSamp)
#print(concatTrueTransp)

#for index, row in concatPreSamp.iterrows():
#  print(row[0])

output=list(reduce(set.intersection, map(set, [dfsTrueTransp[2].gene, dfsTrueTransp[0].gene])))
for i in range(len(output)):
    output[i]=output[i].replace('TRINITY_','')
print(output) ; print(len(output))

