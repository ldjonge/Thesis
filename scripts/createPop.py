#!/usr/bin/env python3
import random
from standardFunctions import *
from numpy.random import choice
from numpy import exp
import sys
from math import log

phenoDict = {"pp":"A", "pq":"A", "pr":"A", "qp":"A", "rp":"A", "qq":"I", "qr":"I", "rq":"I", "rr":"O"}

#Starting parameters
N=500 #Population size, should be possible to keep small
K=200 #Carrying capacity
p=0.33 #Allele frequency
q=0.33 #Allele frequency
r=0.34 #Allele frequency
nGen = 10000 #Number of generations

s = -3 #Selection Pressure
equilA = 0.5 #Equilibrium frequency for andromorphs
equilI = 0.3 #Equilibrium frequency for infuscans

class Male:
    def __init__(self, pAll, mAll):
        self.genotype = pAll+mAll
        self.phenotype = "M"
    def calcFec(self, phenoFreq):
        self.fecundity = 1

class Female:
    def __init__(self,pAll,mAll):
        self.genotype = pAll+mAll
        self.phenotype = phenoDict[self.genotype]

    def calcFec(self, phenoFreq): #As per Le Rouzic 2015
        if self.phenotype == "O":
            self.fecundity = exp(0)
        elif self.phenotype == "I":
            self.fecundity = exp(s*(phenoFreq["I"] - equilI))
        elif self.phenotype == "A":
            self.fecundity = exp(s*(phenoFreq["A"] - equilA))


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
    maleFec = [ind.fecundity for ind in pop[0]]
    maleFecSum = sum(maleFec)
    maleFec = [i*(log(K/maleFecSum)+1) for i in maleFec]
    femaleFec = [ind.fecundity for ind in pop[1]]
    femaleFecSum = sum(femaleFec)

    for i in pop[1]:
        i.fullFec = i.fecundity*(1+log(K/femaleFecSum))
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
    maleFec = [ind.fecundity for ind in pop[0]]
    femaleFec = [ind.fecundity for ind in pop[1]]
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
        fecundity = 1/3/phenoFreq[phenotype]
    except ZeroDivisionError:
        fecundity = N/3
    return fecundity

if __name__ == "__main__":
    malePop = createMalePop(N//2, p,q,r)
    femalePop = createFemalePop(N//2, p,q,r)
    phenFreq = calcPhenoFreq((malePop, femalePop))
    for ind in malePop:
        ind.calcFec(phenFreq)
    for ind in femalePop:
        ind.calcFec(phenFreq)
    pop = (malePop, femalePop)
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
