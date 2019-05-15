import numpy
import pandas
import random

class Gene:
    def __init__(self, var):
        self.gene = []
        if type(var) == int:
            self.geneRange = tuple((0,1) for i in range(var))
        elif type(var) == float:
            try:
                var = int(var)
                self.geneRange = tuple((0,1) for i in range(var))
            except ValueError:
                pass
        elif type(var) == tuple:
            self.geneRange = var
        elif type(var) == list:
            self.geneRange = tuple(var)
        else:
            pass

    def createGene(self):
        for i in self.geneRange:
            self.gene.append(random.uniform(min(i), max(i)))

    def mutate(self, mu):
        for i in range(len(self.gene)):
            if random.random() <= mu:
                self.gene[i] = random.uniform(min(self.geneRange[i]), max(self.geneRange[i]))

    def rangeControl(self):
        for i in range(len(self.gene)):
            if self.gene[i] < min(self.geneRange[i]):
                self.gene[i] = 2*min(self.geneRange[i])-self.gene[i]
            elif self.gene[i] > max(self.geneRange[i]):
                self.gene[i] = 2*max(self.geneRange[i])-self.gene[i]

    def repGene(self, parent1, parent2, id):
        if id == 1:
            self.gene = [0.5*parent1[j]+0.5*parent2[j] for j in range(len(self.geneRange))]
        if id == 2:
            self.gene = [1.5*parent1[j]-0.5*parent2[j] for j in range(len(self.geneRange))]
        if id == 3:
            self.gene = [1.5*parent2[j]-0.5*parent1[j] for j in range(len(self.geneRange))]

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

def crossovers(parents, numVars):
    i=0 #Looping over parents
    offspring = numpy.zeros((1,numVars))
    while i < len(parents):
        parent1 = parents[i].gene
        parent2 = parents[i+1].gene
        offspring = numpy.append(offspring, [[0.5*parent1[j]+0.5*parent2[j] for j in range(numVars)]], axis = 0)
        offspring = numpy.append(offspring, [[1.5*parent1[j]-0.5*parent2[j] for j in range(numVars)]], axis = 0)
        offspring = numpy.append(offspring, [[1.5*parent2[j]-0.5*parent1[j] for j in range(numVars)]], axis = 0)
        i += 2
    offspring = numpy.delete(offspring, 0, axis=0)
    return offspring
def crossover(parents, rangeList):
    i=0
    offspring = []
    while i < len(parents):
        parent1 = parents[i].gene
        parent2 = parents[i+1].gene
        offspring.append(Gene(rangeList))
        offspring[-1].repGene(parent1, parent2, 1)
        offspring.append(Gene(rangeList))
        offspring[-1].repGene(parent1, parent2, 2)
        offspring.append(Gene(rangeList))
        offspring[-1].repGene(parent1, parent2, 3)
        i+= 2
    return offspring

def mutate(offspring, mu):
    mutFree = len(offspring)//4
    for ind in offspring[mutFree:]:
        ind.mutate(mu)
        """
        for i in range(len(ind)):
            if random.random() <= 0.01:
                numpy.put(ind, range(i,i+1), random.random())
        """
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
