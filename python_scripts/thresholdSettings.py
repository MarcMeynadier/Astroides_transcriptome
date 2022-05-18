"""
threshold_settings : Definition of threshold values to process the results of the functionalGenesAnalysis program.

Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                            Analysis computation                              #
#------------------------------------------------------------------------------#



def setThreshold():
    """
    Description
    -----------
    Allows the user to choose the p-value threshold value, set to 0.05 by default.

    Parameters
    ----------
    None

    Returns
    -------
    threshold_pvalue
        float, p-value threshold value.
    """   

    print("\nDefine your threshold value (usually 0.05)\n")
    while True:
        try:
            threshold_pvalue=float(input())
        except ValueError:
            print("\nYou must indicate a valid float ranging from 0 to 1\n")
            continue
        break
    return threshold_pvalue

def filterByCandidate():
    """
    Description
    -----------
    Allows the user to activate the threshold by candidate genes.

    Parameters
    ----------
    None

    Returns
    -------
    flagCandidate
        str, Y or N (yes or no).
    """   

    print("\nDo you want to filter genes by candidate genes ? Y or N\n")
    flagCandidate = input()
    if flagCandidate == 'Y':
        return flagCandidate
    elif flagCandidate == 'N':
        return flagCandidate
    else:
        print('\nYou did not choose a valid answer, no filtration by candidate genes will occur\n')
        flagCandidate = 'N'
        return flagCandidate
