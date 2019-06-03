#!/usr/bin/env python3
import random
from standardFunctions import *
from numpy.random import choice
from numpy import exp
import sys
import math

phenoDict = {"pp":"A", "pq":"A", "pr":"A", "qp":"A", "rp":"A", "qq":"I", "qr":"I", "rq":"I", "rr":"O"}

#Starting parameters
def readParams():
    with open("paramFiles/parameters.ini", "r") as infile:
        paramDict = {}
        for line in infile:
            line = line.strip()
            line = line.split(" ")
            try:
                paramDict[line[0]] = float(line[1])
            except ValueError:
                paramDict[line[0]] = line[1]
    return(paramDict)

"""
N=500 #Population size, should be possible to keep small
K=200 #Carrying capacity
p=0.23 #Allele frequency
q=0.23 #Allele frequency
r=0.54 #Allele frequency
nGen = 100 #Number of generations
s = -3 #Selection Pressure
equilA = 0.5 #Equilibrium frequency for andromorphs
equilI = 0.3 #Equilibrium frequency for infuscans
"""

class Male:
    def __init__(self, pAll, mAll):
        self.genotype = pAll+mAll
        self.phenotype = "M"
    def calcFec(self, popInfo=None):
        if popInfo == None:
            self.fertility = 1
        else:
            self.fertility = popInfo["Mfer"]

    def learning(self, phenoFreq):
        self.aPref = (exp(5*phenoFreq["A"]-2))/(exp(5*phenoFreq["A"]-2)+1)
        self.iPref = (exp(5*phenoFreq["I"]-2))/(exp(5*phenoFreq["I"]-2)+1)
        self.oPref = (exp(5*phenoFreq["O"]-2))/(exp(5*phenoFreq["O"]-2)+1)

    def complexLearning(self, pop):
        if len(pop) >= 50:
            selected = choice(pop, size=50, replace=False)
        else:
            selected = pop
        prefDict = calcPhenoFreqM(selected)
        self.aPref = 0.839*prefDict["M"]+1.872*prefDict["A"]
        self.iPref = 1.898*prefDict["I"]+0.428
        self.oPref = 1.799*prefDict["O"]+0.428

    def genAlgLearning(self, pop, genVars):
        if len(pop) >= 50:
            selected = choice(pop, size=50, replace=False)
        else:
            selected = pop
        prefDict = calcPhenoFreqM(selected)
        self.aPref = genVars[4]*prefDict["M"]+genVars[1]*prefDict["A"]
        self.iPref = genVars[2]*prefDict["I"]+genVars[5]
        self.oPref = genVars[3]*prefDict["O"]+genVars[5]

    def __str__(self):
        return self.phenotype

class Female:
    def __init__(self,pAll,mAll):
        self.genotype = pAll+mAll
        self.phenotype = phenoDict[self.genotype]
        self.taken = 0
        self.mates = []
        self.migrate = 0
        self.migrated = 0

    def __str__(self):
        return self.phenotype

    def calcFec(self, popInfo=None):
        if popInfo == None:
            self.fertility = 1
            self.fecundity = 1
        else:
            if self.phenotype == "O":
                self.fertility = popInfo["Ofer"]
                self.fecundity = popInfo["Ofec"]
            elif self.phenotype == "I":
                self.fertility = popInfo["Ifer"]
                self.fecundity = popInfo["Ifec"]
            elif self.phenotype == "A":
                self.fertility = popInfo["Afer"]
                self.fecundity = popInfo["Afec"]

    def mate(self, male):
        self.taken += 2
        if len(self.mates) < 1:
            self.fecundity *=1.1
        else:
            self.fecundity *=0.9
        self.mates.append(male.genotype)

    def eggLay(self, pop, nEggs):
        nEggs = nEggs * self.fecundity
        nMates = len(self.mates)
        if nMates > 0:
            for male in self.mates[:-1]:
                for i in range(int(0.265/(nMates-1)*nEggs)//8):
                    pop.append(Male(male[0], self.genotype[0]))
                    pop.append(Male(male[0], self.genotype[1]))
                    pop.append(Male(male[1], self.genotype[0]))
                    pop.append(Male(male[1], self.genotype[1]))
                    pop.append(Female(male[0], self.genotype[0]))
                    pop.append(Female(male[0], self.genotype[1]))
                    pop.append(Female(male[1], self.genotype[0]))
                    pop.append(Female(male[1], self.genotype[1]))
            male = self.mates[-1]
            for i in range(int(0.735*nEggs)//8):
                pop.append(Male(male[0], self.genotype[0]))
                pop.append(Male(male[0], self.genotype[1]))
                pop.append(Male(male[1], self.genotype[0]))
                pop.append(Male(male[1], self.genotype[1]))
                pop.append(Female(male[0], self.genotype[0]))
                pop.append(Female(male[0], self.genotype[1]))
                pop.append(Female(male[1], self.genotype[0]))
                pop.append(Female(male[1], self.genotype[1]))
        elif nMates == 1:
            male = self.mates[-1]
            for i in range(nEggs//8):
                pop.append(Male(male[0], self.genotype[0]))
                pop.append(Male(male[0], self.genotype[1]))
                pop.append(Male(male[1], self.genotype[0]))
                pop.append(Male(male[1], self.genotype[1]))
                pop.append(Female(male[0], self.genotype[0]))
                pop.append(Female(male[0], self.genotype[1]))
                pop.append(Female(male[1], self.genotype[0]))
                pop.append(Female(male[1], self.genotype[1]))
        else:
            pass

    def genAlgEggLay(self, pop, nEggs,  genVars):
        nMates = len(self.mates)
        if nMates > 0:
            for male in self.mates[:-1]:
                for i in range(int((1-genVars[12])/(nMates-1)*nEggs)//8):
                    pop.append(Male(male[0], self.genotype[0]))
                    pop.append(Male(male[0], self.genotype[1]))
                    pop.append(Male(male[1], self.genotype[0]))
                    pop.append(Male(male[1], self.genotype[1]))
                    pop.append(Female(male[0], self.genotype[0]))
                    pop.append(Female(male[0], self.genotype[1]))
                    pop.append(Female(male[1], self.genotype[0]))
                    pop.append(Female(male[1], self.genotype[1]))
            male = self.mates[-1]
            for i in range(int(genVars[12]*nEggs)//8):
                pop.append(Male(male[0], self.genotype[0]))
                pop.append(Male(male[0], self.genotype[1]))
                pop.append(Male(male[1], self.genotype[0]))
                pop.append(Male(male[1], self.genotype[1]))
                pop.append(Female(male[0], self.genotype[0]))
                pop.append(Female(male[0], self.genotype[1]))
                pop.append(Female(male[1], self.genotype[0]))
                pop.append(Female(male[1], self.genotype[1]))
        elif nMates == 1:
            male = self.mates[-1]
            for i in range(nEggs//8):
                pop.append(Male(male[0], self.genotype[0]))
                pop.append(Male(male[0], self.genotype[1]))
                pop.append(Male(male[1], self.genotype[0]))
                pop.append(Male(male[1], self.genotype[1]))
                pop.append(Female(male[0], self.genotype[0]))
                pop.append(Female(male[0], self.genotype[1]))
                pop.append(Female(male[1], self.genotype[0]))
                pop.append(Female(male[1], self.genotype[1]))
        else:
            pass

def randomAllele(p,q,r):
    newRand = random.random()
    if newRand <= p:
        allele = "p"
    elif newRand <= p+q:
        allele = "q"
    else:
        allele = "r"
    return allele

def createMalePop(N, p, q, r):
    population = []
    pAll = randomAllele(p,q,r)
    mAll = randomAllele(p,q,r)
    for i in range(N):
        pAll = randomAllele(p,q,r)
        mAll = randomAllele(p,q,r)
        population.append(Male(pAll,mAll))
    return(population)

def createFemalePop(N, p, q, r):
    population = []
    pAll = randomAllele(p,q,r)
    mAll = randomAllele(p,q,r)
    for i in range(N):
        pAll = randomAllele(p,q,r)
        mAll = randomAllele(p,q,r)
        population.append(Female(pAll,mAll))
    return(population)

def reproduce(male, female):
    if random.random() <= 0.5:
        pAll = male.genotype[0]
    else:
        pAll = male.genotype[1]
    if random.random() <= 0.5:
        mAll = female.genotype[0]
    else:
        mAll = female.genotype[1]
    if random.random() <= 0.5:
        return Male(pAll, mAll)
    else:
        return Female(pAll, mAll)

def repVarPop(pop, K):
    newMalePop = []
    newFemalePop = []
    phenoFreq = calcPhenoFreq(pop)
    maleFec = [ind.fertility for ind in pop[0]]
    maleFecSum = sum(maleFec)
    maleFec = [i*(math.log(K/maleFecSum)+1) for i in maleFec]
    femaleFec = [ind.fertility for ind in pop[1]]
    femaleFecSum = sum(femaleFec)

    for i in pop[1]:
        i.fullFec = i.fertility*(1+math.log(K/femaleFecSum))
        while i.fullFec > 1:
            pat = choice(pop[0], 1, maleFec)[0]
            offspring = reproduce(pat, i)
            if type(offspring) == Male:
                newMalePop.append(offspring)
            else:
                newFemalePop.append(offspring)
            i.fullFec -= 1
        if random.random() <= i.fullFec:
            pat = choice(pop[0], 1, maleFec)[0]
            offspring = reproduce(pat, i)
            if type(offspring) == Male:
                newMalePop.append(offspring)
            else:
                newFemalePop.append(offspring)

    newPop = (newMalePop, newFemalePop)
    return(newPop)

def newGen(pop):
    newMalePop = []
    newFemalePop = []
    phenoFreq = calcPhenoFreq(pop)
    maleFec = [ind.fertility for ind in pop[0]]
    femaleFec = [ind.fertility for ind in pop[1]]
    for i in range(N):
        pat = choice(pop[0], 1, maleFec)[0]
        mat = choice(pop[1], 1, femaleFec)[0]
        offspring = reproduce(pat, mat)
        if type(offspring) == Male:
            newMalePop.append(offspring)
        else:
            newFemalePop.append(offspring)
    newPop = (newMalePop, newFemalePop)
    return(newPop)

def calcPhenoFreq(pop):
    freqDict = {"M":0, "A":0, "I":0, "O":0}
    for ind in pop[0]:
        freqDict[ind.phenotype] += 1
    for ind in pop[1]:
        freqDict[ind.phenotype] += 1
    phenoDist = {}
    for phen in "AIO":
        try:
            phenoDist[phen] = freqDict[phen]/sum([freqDict["A"], freqDict["I"], freqDict["O"]])
        except ZeroDivisionError:
            phenoDist[phen] = 0
    return phenoDist

def calcPhenoFreqM(pop):
    freqDict = {"M":0, "A":0, "I":0, "O":0}
    for ind in pop:
        freqDict[ind.phenotype] += 1
    phenoDist = {}
    for phen in "MAIO":
        try:
            phenoDist[phen] = freqDict[phen]/sum([freqDict["A"], freqDict["I"], freqDict["O"], freqDict["M"]])
        except ZeroDivisionError:
            phenoDist[phen] = 0
    return phenoDist

def calcGenoFreq(pop):
    freqDict = {}
    alleleDict = {"p":0, "q":0, "r":0}
    for ind in pop[0]:
        if ind.genotype in freqDict.keys():
            freqDict[ind.genotype] += 1
        else:
            freqDict[ind.genotype] = 1
        for allele in ind.genotype:
            alleleDict[allele] += 1
    for ind in pop[1]:
        if ind.genotype in freqDict.keys():
            freqDict[ind.genotype] += 1
        else:
            freqDict[ind.genotype] = 1
        for allele in ind.genotype:
            alleleDict[allele] += 1
    return (freqDict, alleleDict)

def calcNFDF(phenotype, phenoFreq):
    try:
        fertility = 1/3/phenoFreq[phenotype]
    except ZeroDivisionError:
        fertility = N/3
    return fertility

def startingPop(popInfo):
    N = int(popInfo["N"])
    p = popInfo["p"]
    q = popInfo["q"]
    r = popInfo["r"]
    malePop = createMalePop(N//2, p,q,r)
    femalePop = createFemalePop(N//2, p,q,r)
    phenFreq = calcPhenoFreq((malePop, femalePop))
    #print(phenFreq)
    totalPop = malePop+femalePop
    for ind in malePop:
        ind.calcFec(popInfo)
        ind.complexLearning(totalPop)
    for ind in femalePop:
        ind.calcFec(popInfo)
    pop = [malePop, femalePop]
    return pop

def genAlgStartingPop(N,p,q,r, genVars):
    malePop = createMalePop(N//2, p,q,r)
    femalePop = createFemalePop(N//2, p,q,r)
    phenFreq = calcPhenoFreq((malePop, femalePop))
    #print(phenFreq)
    totalPop = malePop+femalePop
    for ind in malePop:
        ind.calcFec(phenFreq)
        ind.genAlgLearning(totalPop, genVars)
    for ind in femalePop:
        ind.calcFec(phenFreq)
    pop = (malePop, femalePop)
    return pop

if __name__ == "__main__":
    pop = startingPop(N,p,q,r)
    freqTable = [["A", "I", "O", "Males", "Females"]]
    freqTable.append([str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1]))])
    for gen in range(nGen):
        pop = repVarPop(pop, K)
        phenFreq = calcPhenoFreq(pop)
        for ind in pop[0]:
            ind.calcFec(phenFreq)
        for ind in pop[1]:
            ind.calcFec(phenFreq)
        #print([str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1]))])
        freqTable.append([str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1]))])
    #genoFreq = calcGenoFreq(pop)
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'w') as outfile:
            for line in freqTable:
                print("\t".join(line), file=outfile)
    else:
        for line in freqTable:
            print("\t".join(line))
    #print(genoFreq)
