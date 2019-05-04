import pandas
import numpy
from createPop import *
from mating import eggLay, newPopSize, popControl
from genAlg import *
import sys

"""
Starting frequencies can be fixed at 1/3 for these sims
Parameters to vary:
0   successRate: range(0,0.5)
1:4 learning parameters:    aPref, iPref, oPref, mPref, range(0,2)
5   Heterochrome Bonus preference, range(-0.5,0.5)
6   female affected by mating, range(0.5,1)
7   female affected by failed mating, range(0.5,1)
8   male affected by mating, range(0.5,1)
9   male affected by failed mating, range(0.5,1)
10  pref effect failed mating, range(0.5,1)
11  pref effect successful mating, range(1,2)

    Total gene length: 12
"""

def createGenes():
    succesVars = numpy.random.uniform(low=0, high=0.5, size=(48,1))
    learningVars = numpy.random.uniform(low=0, high=2, size=(48,4))
    femBonus = numpy.random.uniform(low=-0.5, high=0.5, size=(48,1))
    matingVars = numpy.random.uniform(low=0.5, high=1, size=(48,5))
    mateSuccEff = numpy.random.uniform(low=1, high=2, size=(48,1))
    varPop = numpy.concatenate((succesVars, learningVars, femBonus, matingVars, mateSuccEff), axis=1)
    return varPop

def matingSearch(pop, newPop, genVars):
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
            while len(hitChance)<200:
                hitChance.append(0.1)
            totalPop = pop[0] + pop[1]
            while len(totalPop)<200:
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
                    if random.random() <= genVars[0]:
                        eggLay(newPop, male, mate, math.ceil(400*mate.fecundity))
                        mate.mate()
                        mate.fecundity *= genVars[6]
                        male.fecundity *= genVars[8]
                        if mate.phenotype == "A":
                            male.aPref *=genVars[11]
                        elif mate.phenotype == "I":
                            male.iPref *=genVars[11]
                        elif mate.phenotype == "O":
                            male.oPref *=genVars[11]
                    else:
                        mate.fecundity *=genVars[7]
                        if mate.fecundity <0.5:
                            if random.random()**2>mate.fecundity:
                                mate.fecundity = 0
                        if mate.fecundity <= 0:
                            pop[1].remove(mate)
                        if mate.phenotype == "A":
                            male.aPref *=genVars[10]
                        elif mate.phenotype == "I":
                            male.iPref *=genVars[10]
                        elif mate.phenotype == "O":
                            male.oPref *=genVars[10]
                elif type(mate) == Male:
                    mate.fecundity *= genVars[9]
                    if mate.fecundity < 0.5:
                        if random.random()**2>mate.fecundity:
                            mate.fecundity = 0
                    if mate.fecundity <= 0:
                        pop[0].remove(mate)
                    male.aPref *=genVars[10]
    else:
        pass

def runSim(genVars):
    pop = genAlgStartingPop(200, 0.18, 0.24, 0.58, genVars)
    for gen in range(50):
        nextGen = []
        for i in range(30):
            matingSearch(pop, nextGen, genVars)
            for fem in pop[1]:
                if fem.taken !=0:
                    fem.taken -=1

        size = newPopSize(400, nextGen, 200)
        pop = popControl(nextGen, size)
        totalPop = pop[0]+pop[1]
        phenFreq = calcPhenoFreq(pop)
        for ind in pop[0]:
            ind.calcFec(phenFreq)
            ind.genAlgLearning(totalPop, genVars)
        for ind in pop[1]:
            ind.calcFec(phenFreq)
    results = list(phenFreq.values())
    results.append(len(totalPop))
    return results

def genAlgGeneration(genVars):
    results = []
    for gene in genVars:
        results.append(runSim(gene))
        #print("Hello World")
    return results

def controlVars(genVars):
    for gene in genVars:
        if gene[0] < 0:
            numpy.put(gene,0, 0-gene[i])
        if gene[0] > 0.5:
            numpy.put(gene,0, 1-gene[i])
        for i in range(1,5):
            if gene[i] < 0:
                numpy.put(gene,i, 0-gene[i])
            if gene[i] > 2:
                numpy.put(gene,i, 4-gene[i])
        if gene[5] < -0.5:
            numpy.put(gene,5, -1-gene[i])
        if gene[5] > 0.5:
            numpy.put(gene,5, 1-gene[i])
        for i in range(6,11):
            if gene[i] < 0.5:
                numpy.put(gene,i, 1-gene[i])
            if gene[i] > 1:
                numpy.put(gene,i, 2-gene[i])
        if gene[11] < 1:
            numpy.put(gene,11, 2-gene[11])
        if gene[11] > 2:
            numpy.put(gene,11, 4-gene[11])
    return genVars

if __name__ == "__main__":
    genePop = createGenes()
    output = []
    for gen in range(20):
        results = genAlgGeneration(genePop)
        fitnessFrame = fitnessSort(results)
        parents = selectMating(genePop, fitnessFrame,32)
        output.append(parents[0])
        output.append(fitnessFrame.iloc[0])
        offspring = crossover(parents, 12)
        offspring = mutate(offspring)
        genePop = controlVars(offspring)
        #print("Wow")
    for line in output:
        if len(sys.argv) > 1:
            with open(sys.argv[1], "w") as outfile:
                print(line, outfile)
