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
                hitChance.append(max(0.5*male.aPref, 0))
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
                if random.random() < male.fecundity:
                    mate = choice(totalPop, p=hitChance)
                else:
                    mate = None
                if type(mate)==Female:
                    contacts += 1
                    if random.random() <= successRate:
                        matings += 1
                        eggLay(newPop, male, mate, math.ceil(nEggs*mate.fecundity))
                        mate.mate()
                        mate.fecundity *= 0.85
                        male.fecundity *= 0.95
                        if mate.phenotype == "A":
                            male.aPref *=1.4
                        elif mate.phenotype == "I":
                            male.iPref *=1.4
                        elif mate.phenotype == "O":
                            male.oPref *=1.4
                    else:
                        mate.fecundity *=0.95
                        if mate.fecundity <0.5:
                            if random.random()**2>mate.fecundity:
                                mate.fecundity = 0
                        if mate.fecundity <= 0:
                            pop[1].remove(mate)
                        if mate.phenotype == "A":
                            male.aPref *=0.8
                        elif mate.phenotype == "I":
                            male.iPref *=0.8
                        elif mate.phenotype == "O":
                            male.oPref *=0.8
                elif type(mate) == Male:
                    MMcontacts += 1
                    mate.fecundity *= 0.95
                    if mate.fecundity < 0.5:
                        if random.random()**2>mate.fecundity:
                            mate.fecundity = 0
                    if mate.fecundity <= 0:
                        pop[0].remove(mate)
                    male.aPref *=0.8
    else:
        pass
    return [matings, contacts, MMcontacts]

def eggLay(pop, male, female, nEggs):
    female.mate()
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
    avgPop = oldPopSize/nEggs
    newPopSize = randomRound(avgPop*(exp(0.5*(K-avgPop)/K)))
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
    newPop = (malePop, femalePop)
    return(newPop)

def recordFec(pop):
    if len(pop[0]) > 0:
        totalMFec = 0
        for m in pop[0]:
            totalMFec += m.fecundity
        avgMFec = totalMFec/len(pop[0])
    else:
        avgMFec = 0
    if len(pop[1]) > 0:
        totalFFec = 0
        for f in pop[1]:
            totalFFec += f.fecundity
        avgFFec = totalFFec/len(pop[1])
    else:
        avgFFec = 0
    return [avgMFec, avgFFec]

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
        for i in range(30):
            results = matingSearch(pop, nEggs, paramDict["successRate"], nextGen, paramDict["K"])
            matings += results[0]
            contacts += results[1]
            MMcontacts += results[2]
            for fem in pop[1]:
                if fem.taken != 0:
                    fem.taken -= 1
        size = newPopSize(nEggs, nextGen, paramDict["K"])
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
            freqTable.append([str(gen), str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1])), str(len(pop[0])+len(pop[1])), str(avgFecs[0]), str(avgFecs[1]), str(prefs[0]), str(prefs[1]), str(prefs[2]), str(matings), str(contacts), str(MMcontacts)])
        elif length == "long":
            freqTable.append([gen, "A", phenFreq["A"]])
            freqTable.append([gen, "I", phenFreq["I"]])
            freqTable.append([gen, "O", phenFreq["O"]])
            freqTable.append([gen, "M", len(pop[0])])
            freqTable.append([gen, "F", len(pop[1])])
            freqTable.append([gen, "T", len(pop[0])+len(pop[1])])
            freqTable.append([gen-1, "MalF", avgFecs[0]])
            freqTable.append([gen-1,"FemF", avgFecs[1]])
            freqTable.append([gen-1, "APref", prefs[0]])
            freqTable.append([gen-1, "IPref", prefs[1]])
            freqTable.append([gen-1, "OPref", prefs[2]])
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
