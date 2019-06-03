import createPop
import pandas
from mating import *
from createPop import *
from standardFunctions import *
from copy import copy
import random

def migrate(pops, female):
    listIndex = None
    for pop in pops:
        if female in pop[1]:
            listIndex = pops.index(pop)
    if listIndex == 0:
        pops[1][1].append(copy(female))
    elif listIndex == len(pops)-1:
        pops[-2][1].append(copy(female))
    else:
        if random.random() <= 0.5:
            pops[listIndex-1][1].append(copy(female))
        else:
            pops[listIndex+1][1].append(copy(female))

def readPopInfo(file="paramFiles/popInfo.csv"):
    with open(file, 'r') as data:
        header = data.readline()
        header = header.strip()
        header = header.split(",")
        pops = []
        for line in data.readlines():
            line = line.strip()
            line = line.split(",")
            popDict={}
            for i in range(len(header)):
                popDict[header[i]] = float(line[i])
            pops.append(popDict)
        print(pops[0].keys())
        return pops

def readMigration(file="paramFiles/dispersalMatrix.csv"):
    data = pandas.read_csv(file, sep=",", header=None)
    data = data.div(data.sum(axis=1), axis=0)
    return(data)

def varMigrate(pops, migrationMatrix, gen, nEggs, length, freqTable):
    nextGen = [[] for i in range(len(pops))]
    for pop in pops:
        migrations = 0
        id = pops.index(pop)
        mRates = migrationMatrix.iloc[id]
        mChance = mRates.cumsum()
        for female in pop[1]:
            migOut = random.random()
            for i in range(len(mChance)):
                if migOut <= mChance[i]:
                    newPop = i
                    break
            if newPop != id:
                migrations += 1
            female.eggLay(nextGen[i], nEggs)
        if length == "long":
            freqTable.append([id+1, gen, "migrations", migrations])
    return nextGen


def runSim(length):
    paramDict = readParams()
    #pops = [startingPop(int(paramDict["N"]), paramDict["p"], paramDict["q"], paramDict["r"]) for i in range(int(paramDict['nPops']))]
    #pops = [startingPop(int(paramDict["N"]), paramDict["p"], paramDict["q"], paramDict["r"]), startingPop(int(paramDict["N"])//2, 1, 0, 0)]
    pops = []
    params = readPopInfo()
    migrationMatrix = readMigration()
    nPop = min(len(params), len(migrationMatrix.index))
    for i in range(nPop):
        pop = params[i]
        pops.append(startingPop(pop))
    nEggs = paramDict["nEggs"]
    if length=="short":
        freqTable = [["Pop", "Gen", "A", "I", "O", "Males", "Females", "Total", "MaleFec", "FemFec", "APref", "IPref", "OPref", "Matings", "Contacts", "MMContacts"]]
    elif length=="long":
        freqTable = [["Pop", "Gen", "Pheno", "Value"]]
    for gen in range(1, int(paramDict["nGen"])+1):
        totalLen=0
        for pop in pops:
            id = pops.index(pop)+1
            nextGen = []
            matings = 0
            contacts = 0
            MMcontacts = 0
            if length == "long":
                prefs = recordPref(pop)
                freqTable.append([id, gen, "preAPref", prefs[0]])
                freqTable.append([id, gen, "preIPref", prefs[1]])
                freqTable.append([id, gen, "preOPref", prefs[2]])
            for i in range(30):
                results = matingSearch(pop, nEggs, paramDict["successRate"], nextGen, int(params[id-1]["K"]))
                matings += results[0]
                contacts += results[1]
                MMcontacts += results[2]
                for fem in pop[1]:
                    if fem.taken != 0:
                        fem.taken -= 1
            avgFecs = recordFec(pop)
            prefs = recordPref(pop)
            if length == "long":
                freqTable.append([id, gen, "Matings", matings])
                freqTable.append([id, gen, "Contacts", contacts])
                freqTable.append([id, gen, "MMContacts", MMcontacts])
                freqTable.append([id, gen, "MalF", avgFecs[0]])
                freqTable.append([id, gen,"FemF", avgFecs[1]])
                freqTable.append([id, gen, "APref", prefs[0]])
                freqTable.append([id, gen, "IPref", prefs[1]])
                freqTable.append([id, gen, "OPref", prefs[2]])
        eggs = varMigrate(pops, migrationMatrix, gen, nEggs, length, freqTable)
        newPops = []
        for pop in eggs:
            id = eggs.index(pop)+1
            size = newPopSize(nEggs, pop, int(params[id-1]["K"]))
            newPop = popControl(pop, size)
            totalPop = newPop[0]+newPop[1]
            totalLen += len(totalPop)
            phenFreq = calcPhenoFreq(newPop)
            for ind in newPop[0]:
                ind.calcFec(params[id-1])
                ind.complexLearning(totalPop)
            for ind in newPop[1]:
                ind.calcFec(params[id-1])
            #if length == "short":
            #    freqTable.append([str(id), str(gen), str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1])), str(len(pop[0])+len(pop[1])), str(avgFecs[0]), str(avgFecs[1]), str(prefs[0]), str(prefs[1]), str(prefs[2]), str(matings), str(contacts), str(MMcontacts)])
            if length == "long":
                freqTable.append([id, gen, "A", phenFreq["A"]])
                freqTable.append([id, gen, "I", phenFreq["I"]])
                freqTable.append([id, gen, "O", phenFreq["O"]])
                freqTable.append([id, gen, "M", len(newPop[0])])
                freqTable.append([id, gen, "F", len(newPop[1])])
                freqTable.append([id, gen, "T", len(newPop[0])+len(newPop[1])])
            else:
                print("failure")
            newPops.append(newPop)
            if len(pops[id-1][0])+len(pops[id-1][1]) == 0 and len(totalPop) > 0:
                print("New population formed in pop {}, generation {}".format(str(id), str(gen)))
            if len(pops[id-1][0])+len(pops[id-1][1]) > 0 and len(totalPop) == 0:
                print("Population {} extinct in generation {}".format(str(id), str(gen)))
        if (gen)%10 == 0:
            print("Generation {} complete, population size {}".format(str(gen), str(totalLen)))
        pops = newPops
    return freqTable

if __name__ == "__main__":
    freqTable = runSim("long")
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'w') as outfile:
            for line in freqTable:
                print("\t".join(line), file=outfile)
    else:
        for line in freqTable:
            print("\t".join(line))
