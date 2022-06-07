"""
getCandidates : From the annotation file built with hmmsearch using the transcriptome, retrieves the genes corresponding 
to the candidate genes from a set of keywords specific to the biological problem of the project.

Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                       Importation of external libraries                      #
#------------------------------------------------------------------------------#


import pandas as pd


#------------------------------------------------------------------------------#
#                              Files management                                #
#------------------------------------------------------------------------------#


def getAnnotation(): 
    """
    Description
    -----------
    Retrieves the output file from hmmsearch, parse it to retrieve the Pfam annotation and the Pfam code of the genes. 
    Retrieves the dataframe from getProteinSequences(), then the output file from hmmsearch in order to parse it to get 
    the Pfam annotation as well as the Pfam code of the genes. The two data frames are then merged together and returned 
    as a common data frame. The merge is done so that even a non-coding gene ID is returned. 

    Parameters
    ----------
    None

    Returns
    -------
    mergeDf
       Pandas dataframe, contains successively the gene IDs, the Pfam annotations, the Pfam codes and the protein sequences. 
    """

    curedFile = []
    with open('../../data/net/8_functionalAnnotation/hmmsearchOutput.out') as f:
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


#------------------------------------------------------------------------------#
#                            Analysis computation                              #
#------------------------------------------------------------------------------#


def getCandidateGenes():
    """
    Description
    -----------
    From the annotation file and a dictionary containing keywords specific to proteins involved in the calcification process. 
    Duplicate genes are removed and the cleaned up candidate gene file is then saved. 

    Parameters
    ----------
    None

    Returns
    -------
    None 
    """

    annot = getAnnotation()
    candidate=['lectin','sushi','galaxin','collagen','adhesin','cadherin','actin','HCO3','anhydrase','V_ATPase','calmodulin']
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
    
