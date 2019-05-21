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

def runSim(length):
    paramDict = readParams()
    #pops = [startingPop(int(paramDict["N"]), paramDict["p"], paramDict["q"], paramDict["r"]) for i in range(int(paramDict['nPops']))]
    pops = [startingPop(int(paramDict["N"]), paramDict["p"], paramDict["q"], paramDict["r"]), startingPop(int(paramDict["N"])//2, 1, 0, 0)]
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
                results = matingSearch(pop, nEggs, paramDict["successRate"], nextGen, paramDict["K"])
                matings += results[0]
                contacts += results[1]
                MMcontacts += results[2]
                for fem in pop[1]:
                    if fem.taken != 0:
                        fem.taken -= 1
            if length == "long":
                freqTable.append([id, gen, "Matings", matings])
                freqTable.append([id, gen, "Contacts", contacts])
                freqTable.append([id, gen, "MMContacts", MMcontacts])
        for pop in pops:
            for fem in pop[1]:
                if random.random() <= float(paramDict["m"]) and fem.migrate == 0:
                    fem.migrate += 1
                    migrate(pops, fem)
                    fem.migrated += 1
                if fem.migrate > 1:
                    print("Wow")
        newPops = []
        for pop in pops:
            nextGen = []
            id = pops.index(pop)+1
            for fem in pop[1]:
                if fem.migrated == 0:
                    fem.eggLay(nextGen, nEggs)
            size = newPopSize(nEggs, nextGen, int(paramDict["K"]))
            avgFecs = recordFec(pop)
            prefs = recordPref(pop)

            newPop = popControl(nextGen, size)
            totalPop = newPop[0]+newPop[1]
            totalLen += len(totalPop)
            phenFreq = calcPhenoFreq(newPop)
            for ind in newPop[0]:
                ind.calcFec(phenFreq)
                ind.complexLearning(totalPop)
            for ind in newPop[1]:
                ind.calcFec(phenFreq)
            if length == "short":
                freqTable.append([str(id), str(gen), str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1])), str(len(pop[0])+len(pop[1])), str(avgFecs[0]), str(avgFecs[1]), str(prefs[0]), str(prefs[1]), str(prefs[2]), str(matings), str(contacts), str(MMcontacts)])
            elif length == "long":
                freqTable.append([id, gen, "A", phenFreq["A"]])
                freqTable.append([id, gen, "I", phenFreq["I"]])
                freqTable.append([id, gen, "O", phenFreq["O"]])
                freqTable.append([id, gen, "M", len(newPop[0])])
                freqTable.append([id, gen, "F", len(newPop[1])])
                freqTable.append([id, gen, "T", len(newPop[0])+len(newPop[1])])
                freqTable.append([id, gen, "MalF", avgFecs[0]])
                freqTable.append([id, gen,"FemF", avgFecs[1]])
                freqTable.append([id, gen, "APref", prefs[0]])
                freqTable.append([id, gen, "IPref", prefs[1]])
                freqTable.append([id, gen, "OPref", prefs[2]])
            else:
                print("failure")
            newPops.append(newPop)
        if (gen)%10 == 0:
            print("Generation {} complete, population size {}".format(str(gen), str(totalLen)))
        pops = newPops
    return freqTable

if __name__ == "__main__":
    freqTable = runSim("short")
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'w') as outfile:
            for line in freqTable:
                print("\t".join(line), file=outfile)
    else:
        for line in freqTable:
            print("\t".join(line))
