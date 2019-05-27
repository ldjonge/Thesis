import random
from standardFunctions import *
from numpy.random import choice
from numpy import exp
import sys
import math
from createPop import *

def matingSearch(pop, nEggs, successRate, newPop, K):
    matings = 0
    contacts = 0
    MMcontacts = 0
    if len(pop[1])!=0:
        for male in pop[0]:
            hitChance = []
            for ind in pop[0]:
                hitChance.append(max(0.0796*male.aPref, 0))
            for fem in pop[1]:
                if fem.taken != 0:
                    hitChance.append(0)
                elif fem.phenotype == "A":
                    hitChance.append(max(male.aPref, 0))
                elif fem.phenotype == "I":
                    hitChance.append(max(male.iPref, 0))
                elif fem.phenotype == "O":
                    hitChance.append(max(male.oPref, 0))
            while len(hitChance)<K:
                hitChance.append(0.1)
            totalPop = pop[0] + pop[1]
            while len(totalPop)<K:
                totalPop.append(None)
            totalW = sum(hitChance)
            if totalW != 0:
                for i in range(len(hitChance)):
                    hitChance[i] = hitChance[i]/totalW
                if random.random() < male.fertility:
                    mate = choice(totalPop, p=hitChance)
                else:
                    mate = None
                if type(mate)==Female:
                    contacts += 1
                    if random.random() <= successRate*mate.fertility:
                        matings += 1
                        mate.mate(male)
                        mate.fertility *= 0.896
                        male.fertility *= 0.963
                        if mate.phenotype == "A":
                            male.aPref *=1.179
                        elif mate.phenotype == "I":
                            male.iPref *=1.179
                        elif mate.phenotype == "O":
                            male.oPref *=1.179
                    else:
                        mate.fertility *=0.977
                        if mate.fertility <0.5:
                            if random.random()**2>mate.fertility:
                                mate.fertility = 0
                        #if mate.fertility <= 0:
                        #    pop[1].remove(mate)
                        if mate.phenotype == "A":
                            male.aPref *=0.921
                        elif mate.phenotype == "I":
                            male.iPref *=0.921
                        elif mate.phenotype == "O":
                            male.oPref *=0.921
                elif type(mate) == Male:
                    MMcontacts += 1
                    mate.fertility *= 0.821
                    if mate.fertility < 0.5:
                        if random.random()**2>mate.fertility:
                            mate.fertility = 0
                    #if mate.fertility <= 0:
                    #    pop[0].remove(mate)
                    male.aPref *=0.921
    else:
        pass
    return [matings, contacts, MMcontacts]


def eggLay(pop, male, female, nEggs):
    female.mate(male)
    for i in range(nEggs//8):
        pop.append(Male(male.genotype[0], female.genotype[0]))
        pop.append(Male(male.genotype[0], female.genotype[1]))
        pop.append(Male(male.genotype[1], female.genotype[0]))
        pop.append(Male(male.genotype[1], female.genotype[1]))
        pop.append(Female(male.genotype[0], female.genotype[0]))
        pop.append(Female(male.genotype[0], female.genotype[1]))
        pop.append(Female(male.genotype[1], female.genotype[0]))
        pop.append(Female(male.genotype[1], female.genotype[1]))

def newPopSize(nEggs, pop, K):
    oldPopSize = len(pop)
    avgPop = 2*oldPopSize/nEggs
    newPopSize = randomRound(avgPop*(exp(0.5*(K-avgPop)/K)))
    if avgPop > K:
        newPopSize = max(K, newPopSize)
    return(newPopSize)

def popControl(pop, size):
    newPop = list(choice(pop, size=size, replace=False))
    malePop = []
    femalePop = []
    for ind in newPop:
        if type(ind)==Male:
            malePop.append(ind)
        elif type(ind)==Female:
            femalePop.append(ind)
        else:
            print("error, individual of type: ", type(ind))
    newPop = [malePop, femalePop]
    return(newPop)

def recordFec(pop):
    if len(pop[0]) > 0:
        totalMFec = 0
        for m in pop[0]:
            totalMFec += m.fertility
        avgMFec = totalMFec/len(pop[0])
    else:
        avgMFec = 0
    if len(pop[1]) > 0:
        totalFFec = 0
        totalAFec = 0
        totalIFec = 0
        totalOFec = 0
        totalA = 0
        totalI = 0
        totalO = 0
        for f in pop[1]:
            totalFFec += f.fertility
            if f.phenotype == "A":
                totalAFec += f.fertility
                totalA += 1
            elif f.phenotype == "I":
                totalIFec += f.fertility
                totalI += 1
            elif f.phenotype == "O":
                totalOFec += f.fertility
                totalO += 1
        avgFFec = totalFFec/len(pop[1])
        try:
            avgAFec = totalAFec/totalA
        except ZeroDivisionError:
            avgAFec = 0
        try:
            avgIFec = totalIFec/totalI
        except ZeroDivisionError:
            avgIFec = 0
        try:
            avgOFec = totalOFec/totalO
        except ZeroDivisionError:
            avgOFec = 0
    else:
        avgFFec = 0
        avgAFec = 0
        avgIFec = 0
        avgOFec = 0
    return [avgMFec, avgFFec, avgAFec, avgIFec, avgOFec]

def recordPref(pop):
    totalAPref = 0
    totalIPref = 0
    totalOPref = 0
    for m in pop[0]:
        totalAPref += m.aPref
        totalIPref += m.iPref
        totalOPref += m.oPref
    prefSum = totalAPref + totalIPref + totalOPref
    if prefSum > 0:
        avgAPref = totalAPref/prefSum
        avgIPref = totalIPref/prefSum
        avgOPref = totalOPref/prefSum
        return [avgAPref, avgIPref, avgOPref]
    else:
        return [0,0,0]


def runSim(length):
    paramDict = readParams()
    pop = startingPop(int(paramDict["N"]), paramDict["p"], paramDict["q"], paramDict["r"])
    nEggs = paramDict["nEggs"]
    if length=="short":
        freqTable = [["Gen", "A", "I", "O", "Males", "Females", "Total", "MaleFec", "FemFec", "APref", "IPref", "OPref", "Matings", "Contacts", "MMContacts"]]
    elif length=="long":
        freqTable = [["Gen", "Pheno", "Value"]]
    for gen in range(1, int(paramDict["nGen"])+1):
        nextGen = []
        matings = 0
        contacts = 0
        MMcontacts = 0
        if length == "long":
            prefs = recordPref(pop)
            freqTable.append([gen, "preAPref", prefs[0]])
            freqTable.append([gen, "preIPref", prefs[1]])
            freqTable.append([gen, "preOPref", prefs[2]])
        for i in range(30):
            results = matingSearch(pop, nEggs, paramDict["successRate"], nextGen, paramDict["K"])
            matings += results[0]
            contacts += results[1]
            MMcontacts += results[2]
            for fem in pop[1]:
                if fem.taken != 0:
                    fem.taken -= 1
        for fem in pop[1]:
            fem.eggLay(nextGen, nEggs)
        size = newPopSize(nEggs, nextGen, int(paramDict["K"]))
        avgFecs = recordFec(pop)
        prefs = recordPref(pop)

        pop = popControl(nextGen, size)
        totalPop = pop[0]+pop[1]
        phenFreq = calcPhenoFreq(pop)
        for ind in pop[0]:
            ind.calcFec(phenFreq)
            ind.complexLearning(totalPop)
        for ind in pop[1]:
            ind.calcFec(phenFreq)
        #print([str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1]))])
        if length == "short":
            freqTable.append([str(gen), str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1])), str(len(pop[0])+len(pop[1])), str(avgFecs[0]), str(avgFecs[1]), str(avgFecs[2]), str(avgFecs[3]), str(avgFecs[4]), str(prefs[0]), str(prefs[1]), str(prefs[2]), str(matings), str(contacts), str(MMcontacts)])
        elif length == "long":
            freqTable.append([gen, "A", phenFreq["A"]])
            freqTable.append([gen, "I", phenFreq["I"]])
            freqTable.append([gen, "O", phenFreq["O"]])
            freqTable.append([gen, "M", len(pop[0])])
            freqTable.append([gen, "F", len(pop[1])])
            freqTable.append([gen, "T", len(pop[0])+len(pop[1])])
            freqTable.append([gen, "MalF", avgFecs[0]])
            freqTable.append([gen,"FemF", avgFecs[1]])
            freqTable.append([gen,"AF", avgFecs[2]])
            freqTable.append([gen,"IF", avgFecs[3]])
            freqTable.append([gen,"OF", avgFecs[4]])
            freqTable.append([gen, "APref", prefs[0]])
            freqTable.append([gen, "IPref", prefs[1]])
            freqTable.append([gen, "OPref", prefs[2]])
            freqTable.append([gen, "Matings", matings])
            freqTable.append([gen, "Contacts", contacts])
            freqTable.append([gen, "MMContacts", MMcontacts])

        else:
            print("failure")
        if (gen)%10 == 0:
            print("Generation {} complete, population size {}".format(str(gen), str(len(pop[0])+len(pop[1]))))
    return freqTable

def main():
    freqTable = runSim("long")
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'w') as outfile:
            for line in freqTable:
                print("\t".join(line), file=outfile)
    else:
        for line in freqTable:
            print("\t".join(line))

if __name__ == "__main__":
    main()
