import numpy
import pandas
import random

numVars = 3
popSize = 30
fullPop = (popSize, numVars)
newPop = numpy.random.uniform(low=0, high=1, size=fullPop)

def controlTotal(pop, numVars):
    for ind in pop:
        for i in range(numVars):
            if ind[i] < 0:
                numpy.put(ind,range(i,i+1), 0-ind[i])
            if ind[i] > 1:
                numpy.put(ind,range(i,i+1), 2-ind[i])
        total = ind.sum()
        for i in range(numVars):
            numpy.put(ind, range(i,i+1), ind[i]/total)
    return(pop)
newPop = controlTotal(newPop, numVars)

nGen = 1000
nMating = popSize//3*2

def fitnessSort(pop):
    indx = []
    freqFit = []

    for i in range(len(pop)):
        ind = pop[i]
        indx.append(i)
        freqFit.append(ind[0]*ind[1]*ind[2])
    fitnessFrame = pandas.DataFrame({"index": indx, "Freqs": freqFit})
    fitnessFrame.sort_values(by="Freqs", ascending=False, inplace=True)
    return fitnessFrame

def selectMating(newPop, fitnessFrame, nMating):
    parentIndex = fitnessFrame.iloc[0:nMating,0]
    parents = []
    for i in parentIndex:
        parents.append(newPop[i])
    return(parents)

def crossover(parents, numVars):
    i=0 #Looping over parents
    offspring = numpy.zeros((1,numVars))
    while i < len(parents):
        parent1 = parents[i]
        parent2 = parents[i+1]
        offspring = numpy.append(offspring, [[0.5*parent1[j]+0.5*parent2[j] for j in range(numVars)]], axis = 0)
        offspring = numpy.append(offspring, [[1.5*parent1[j]-0.5*parent2[j] for j in range(numVars)]], axis = 0)
        offspring = numpy.append(offspring, [[1.5*parent2[j]-0.5*parent1[j] for j in range(numVars)]], axis = 0)
        i += 2
    offspring = numpy.delete(offspring, 0, axis=0)
    return offspring

def mutate(offspring):
    mutFree = len(offspring)//4
    for ind in offspring[mutFree:]:
        for i in range(len(ind)):
            if random.random() <= 0.01:
                numpy.put(ind, range(i,i+1), random.random())
    return offspring

if __name__ == "__main__":
    for gen in range(nGen):
        fitnessFrame = fitnessSort(newPop)
        parents = selectMating(newPop, fitnessFrame, nMating)
        offspring = crossover(parents, numVars)
        offspring = mutate(offspring)
        newPop = controlTotal(offspring, numVars)

    print(newPop)
    print(fitnessSort(newPop))
