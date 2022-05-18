"""
averageLength_SD_bioanalyser : From the results provided by a bioanalyzer, provide the average length of the RNA fragments 
and the standard deviation of this length. These two metrics are necessary to run the Kallisto tool in single end mode.  

Marc Meynadier
"""


#------------------------------------------------------------------------------#
#                            Analysis computation                              #
#------------------------------------------------------------------------------#


from numpy import mean,std


#------------------------------------------------------------------------------#
#                            Analysis computation                              #
#------------------------------------------------------------------------------#


def calcMetrics():
    """
    Description
    -----------
    From the information provided by the user, calculates a mean and the standard deviation of a set of values. 
    Initially designed to calculate the average length and standard deviation of RNA fragments analyzed by a bioanalyzer, 
    the user can choose to provide the direct length of the fragment or the indirect length (beginning of the fragment and 
    end of the fragment, depending on the bioanalyzer used).

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    print("How many graphs samples are available in your dataset ?")
    graphNumber = int(input())
    print("Fragment length directly or indirectly ? (y/n)")
    typeResult = input()
    fragmentsLengthList = []
    if typeResult == 'n': 
        for i in range(graphNumber):
            print("Beginning of fragment number",i+1,":") ; bFragment = int(input())
            print("End of fragment",i+1,":"); eFragment = int(input())
            fragmentLength = eFragment - bFragment ; fragmentsLengthList.append(fragmentLength)
    elif typeResult == 'y':
        for i in range(graphNumber):
            print("Length of fragment number",i+1,":") ; fragmentLength = int(input())
            fragmentsLengthList.append(fragmentLength)
    averageFragmentLength = mean(fragmentsLengthList)
    sdFragmentLength = std(fragmentsLengthList)
    print("Average fragments length = ",round(averageFragmentLength))
    print("Standard deviation of fragments length = ",round(sdFragmentLength))

calcMetrics()
