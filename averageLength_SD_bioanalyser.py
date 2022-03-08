# Automatic calculation average length fragment and standard deviation of length fragment
# Marc Meynadier

from numpy import mean,std

def calcMetrics():
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
