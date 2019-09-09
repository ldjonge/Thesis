import pandas as pd

"""
Function to read in parameter values from external file, allowing for a potential alternative file
The function returns a dictionary with the values for each parameter
"""
def readParams(file="paramFiles/parameters.ini"):
    with open(file, "r") as infile:
        paramDict = {}
        for line in infile:
            line = line.strip()
            line = line.split(" ") # ini file should be space separated
            """
            Numeric parameters should be converted to floats, other parameters should stay strings.
            In case integers are needed for a specific parameter this conversion can be performed later
            """
            try:
                """
                Numeric parameters should be recorded as such (floats), other parameters should stay strings.
                In case integers are needed for a specific parameter this conversion will be performed later
                """
                paramDict[line[0]] = float(line[1])
            except ValueError:
                paramDict[line[0]] = line[1]
    return(paramDict)

"""
All population specific parameter values such as survival rate and fecundity per sex/morph are specified in a csv file.
This file also includes the starting status of each population, including population size and allele frequency
"""
def readPopInfo(file="paramFiles/popInfo.csv"):
    with open(file, 'r') as data:
        header = data.readline()
        header = header.strip()
        header = header.split(",")
        """
        Using the header line of the file makes it so the parameters can be entered in any order,
        and possible missing parameters do not cause any big issues
        """
        pops = []
        for line in data.readlines():
            line = line.strip()
            line = line.split(",")
            popDict={}
            for i in range(len(header)):
                popDict[header[i]] = float(line[i])
            pops.append(popDict)
        """
        By giving each population its own dictionary of parameter values, it is very easy to index them
        """
        return pops

"""
The migration probability between each population, including the probability of staying in the same population, is specified in a csv file.
These probabilities are transformed to make sure they add up to 1, to make sure no errors occur.
This could be adjusted to only accept percentages and proper probabilities so typos are filtered out
"""
def readMigration(file="paramFiles/dispersalMatrix.csv"):
    data = pd.read_csv(file, sep=",", header=None)
    data = data.div(data.sum(axis=1), axis=0) # Dividing by the sum of the row normalises probabilities
    return(data)
