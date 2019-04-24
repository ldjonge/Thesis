import random
from standardFunctions import *
from numpy.random import choice
from numpy import exp
import sys
import math
from createPop import *

nEggs = 208

def matingSearch(pop, nEggs, successRate, newPop, K):
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
                hitChance.append(0.2)
            totalPop = pop[0] + pop[1]
            while len(totalPop)<K:
                totalPop.append(None)
            totalW = sum(hitChance)
            if totalW != 0:
                for i in range(len(hitChance)):
                    hitChance[i] = hitChance[i]/totalW
                mate = choice(totalPop, p=hitChance)
                if type(mate)==Female:
                    if random.random() <= successRate*mate.fecundity:
                        eggLay(newPop, male, mate, math.ceil(nEggs*mate.fecundity))
                        mate.mate()
                    else:
                        mate.fecundity *=0.8
                        #if mate.fecundity <0.5:
                        #    if random.random()**2>mate.fecundity:
                        #        mate.fecundity = 0
                        if mate.fecundity <= 0:
                            mate.fecundity = 0
    else:
        pass

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
    newPopSize = randomRound(avgPop*(exp((K-avgPop)/K)))
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



def runSim(length):
    paramDict = readParams()
    pop = startingPop(int(paramDict["N"]), paramDict["p"], paramDict["q"], paramDict["r"])
    if length=="short":
        freqTable = [["Gen", "A", "I", "O", "Males", "Females"]]
    elif length=="long":
        freqTable = [["Gen", "Pheno", "Value"]]
    for gen in range(1, nGen+1):
        nextGen = []
        for i in range(20):
            matingSearch(pop, nEggs, paramDict["successRate"], nextGen, paramDict["K"])
            for fem in pop[1]:
                if fem.taken != 0:
                    fem.taken -= 1
        size = newPopSize(nEggs, nextGen, paramDict["K"])
        """
        Printing fecundities for sake of testing

        try:
            testMale = pop[0][0]
            print(testMale.aPref, testMale.iPref, testMale.oPref)
        except IndexError:
            print("No Males")

        for fem in pop[1]:
            print(fem.phenotype, ":", fem.fecundity, "\t", fem.taken)
        """
        pop = popControl(nextGen, size)
        phenFreq = calcPhenoFreq(pop)
        for ind in pop[0]:
            ind.calcFec(phenFreq)
            ind.learning(phenFreq)
        for ind in pop[1]:
            ind.calcFec(phenFreq)
        #print([str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1]))])
        if length == "short":
            freqTable.append([str(gen), str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1]))])
        elif length == "long":
            freqTable.append([gen, "A", phenFreq["A"]])
            freqTable.append([gen, "I", phenFreq["I"]])
            freqTable.append([gen, "O", phenFreq["O"]])
            freqTable.append([gen, "M", len(pop[0])])
            freqTable.append([gen, "F", len(pop[1])])

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
