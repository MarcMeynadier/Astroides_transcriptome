import pandas as pd

def getAnnotation(): 
    curedFile = []
    with open('../../data/net/8_functionnalAnnotation/hmmsearchOutput.out') as f:
        contents = f.readlines()
    for i in contents:
        if 'TRINITY' in i:
            curedFile.append(i)
    geneNames=[] ; pfamAnnot= [] ; pfamCode = []
    separator = "PF"
    for i in range(len(curedFile)):
        geneNames.append(curedFile[i].split(" - ",1)[0])
        splitDash = curedFile[i].split(" - ",1)[1]
        splitAnnot = splitDash.split(separator,1)[0] ; splitAnnot = splitAnnot.strip()
        splitCode = splitDash.split(separator,1)[1] ; splitCode = splitCode.split(".",1)[0] ; splitCode = separator + splitCode
        pfamAnnot.append(splitAnnot) ; pfamCode.append(splitCode)
    dic={'genes':geneNames,'pfam_annotation':pfamAnnot,'pfam_code':pfamCode}
    mergeDf=pd.DataFrame(dic) 
    mergeDf = mergeDf.replace(to_replace ='(_i).*', value = '', regex = True)
    return mergeDf

def getCandidateGenes():
    annot = getAnnotation()
    candidate=['ectin','sushi','galaxin','collagen','adhesin','cadherin','actin','HCO3','anhydrase','V_ATPase','calmodulin']
    annotFilter = annot[annot.stack().str.contains('|'.join(candidate)).any(level=0)]
    annotFilter = annotFilter.drop_duplicates(subset='genes', keep='first', inplace=False, ignore_index=False)
    pfam = annotFilter['pfam_annotation'].tolist()
    pfamFilter = []
    for i in pfam:
        newName = i.split(' ',1)[1]
        pfamFilter.append(newName)
    annotFilter['pfam_annotation'] = pfamFilter
    annotFilter.to_csv('candidateGenes.csv',index=False,encoding='utf-8')
    return
    