import pandas
import numpy
from createPop import *
from mating import eggLay, newPopSize, popControl
from genAlg import *
import sys

"""
Starting frequencies can be fixed at 1/3 for these sims
Parameters to vary:
0   successRate: range(0,0.3)
1:4 learning parameters:    aPref, iPref, oPref, mPref, range(0,2)
5   Heterochrome Bonus preference, range(-0.5,0.5)
6   female affected by mating, range(0.8,1)
7   female affected by failed mating, range(0.8,1)
8   male affected by mating, range(0.8,1)
9   male affected by failed mating, range(0.8,1)
10  pref effect failed mating, range(0.5,1)
11  pref effect successful mating, range(1,2)
12  sperm competition, range(0.6,1)
13  male recognition, range(0,1)

    Total gene length: 14
"""

geneParamVars = ((0,0.3), (0,2), (0,2),(0,2),(0,2),(-0.5,0.5), (0.8,1),(0.8,1),(0.8,1),(0.8,1),(0.5,1), (1,2),(0.6,1),(0,1))

def createGenes():
    succesVars = numpy.random.uniform(low=0, high=0.3, size=(48,1))
    learningVars = numpy.random.uniform(low=0, high=2, size=(48,4))
    femBonus = numpy.random.uniform(low=-0.5, high=0.5, size=(48,1))
    matingVars = numpy.random.uniform(low=0.8, high=1, size=(48,4))
    mateFailEff = numpy.random.uniform(low=0.5, high=1, size=(48,1))
    mateSuccEff = numpy.random.uniform(low=1, high=2, size=(48,1))
    spermComp = numpy.random.uniform(low=0.6, high=1, size=(48,1))
    maleRec = numpy.random.uniform(low=0, high=1, size=(48,1))
    varPop = numpy.concatenate((succesVars, learningVars, femBonus, matingVars, mateFailEff, mateSuccEff, spermComp, maleRec), axis=1)
    return varPop

def matingSearch(pop, newPop, genVars):
    if len(pop[1])!=0:
        for male in pop[0]:
            hitChance = []
            for ind in pop[0]:
                hitChance.append(max(genVars[13]*male.aPref, 0))
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
                        #eggLay(newPop, male, mate, math.ceil(400))
                        mate.mate(male)
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
    for gen in range(100):
        nextGen = []
        for i in range(30):
            matingSearch(pop, nextGen, genVars)
            for fem in pop[1]:
                if fem.taken !=0:
                    fem.taken -=1
        for fem in pop[1]:
            fem.genAlgEggLay(nextGen, 400, genVars)
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
    return results

def controlVars(genVars):
    for gene in genVars:
        if gene[0] < 0:
            numpy.put(gene,0, 0-gene[0])
        if gene[0] > 0.3:
            numpy.put(gene,0, 0.6-gene[0])
        for i in range(1,5):
            if gene[i] < 0:
                numpy.put(gene,i, 0-gene[i])
            if gene[i] > 2:
                numpy.put(gene,i, 4-gene[i])
        if gene[5] < -0.5:
            numpy.put(gene,5, -1-gene[5])
        if gene[5] > 0.5:
            numpy.put(gene,5, 1-gene[5])
        for i in range(6,10):
            if gene[i] < 0.8:
                numpy.put(gene,i, 1.6-gene[i])
            if gene[i] > 1:
                numpy.put(gene,i, 2-gene[i])
        if gene[10] < 0.5:
            numpy.put(gene,10, 1-gene[10])
        if gene[10] > 1:
            numpy.put(gene,10, 2-gene[10])
        if gene[11] < 1:
            numpy.put(gene,11, 2-gene[11])
        if gene[11] > 2:
            numpy.put(gene,11, 4-gene[11])
        if gene[12] < 0.6:
            numpy.put(gene,12, 1.2-gene[12])
        if gene[12] > 1:
            numpy.put(gene,12, 2-gene[12])
        if gene[13] < 0:
            numpy.put(gene,13, 0-gene[13])
        if gene[13] > 1:
            numpy.put(gene,13,2-gene[13])
    return genVars

if __name__ == "__main__":
    genePop = createGenes()
    output = []
    for gen in range(20):
        results = genAlgGeneration(genePop)
        fitnessFrame = fitnessSort(results)
        parents = selectMating(genePop, fitnessFrame,32)
        for ind in range(5):
            output.append(parents[ind])
            output.append(results[fitnessFrame.iloc[ind,0]])
        offspring = crossover(parents, 14)
        #offspring = mutate(offspring)
        genePop = controlVars(offspring)
        output.append("Run {}".format(str(gen+1)))
        print("Run {} complete".format(str(gen+1)))
    results = genAlgGeneration(genePop)
    fitnessFrame = fitnessSort(results)
    parents = selectMating(genePop, fitnessFrame,32)
    for ind in range(5):
        output.append(parents[ind])
        output.append(results[fitnessFrame.iloc[ind,0]])
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as outfile:
            for line in output:
                print(line, file=outfile)
    else:
        for line in output:
            print(line)
