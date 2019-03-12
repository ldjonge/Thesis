#!/usr/bin/env python3
import random
from standardFunctions import *
from numpy.random import choice
from numpy import exp
import sys

phenoDict = {"pp":"A", "pq":"A", "pr":"A", "qp":"A", "rp":"A", "qq":"I", "qr":"I", "rq":"I", "rr":"O"}

#Starting parameters
N=100 #Population size, should be possible to keep small
p=2/6 #Allele frequency
q=2/6 #Allele frequency
r=1/3 #Allele frequency
nGen = 100 #Number of generations

s = -2 #Selection Pressure
equilA = 0.6 #Equilibrium frequency for andromorphs
equilI = 0.15 #Equilibrium frequency for infuscans

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
    '''
    for phen in ["M", "A", "I", "O"]:
        print("{}:\t{}\t{}%".format(phen, freqDict[phen], 100*freqDict[phen]/sum(freqDict.values())))
    '''
    phenoDist = {}
    for phen in "AIO":
        phenoDist[phen] = freqDict[phen]/sum([freqDict["A"], freqDict["I"], freqDict["O"]])
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
    '''
    for gen in freqDict.keys():
        print("{}:\t{}\t{}%".format(gen, freqDict[gen], 200*freqDict[gen]/sum(freqDict.values())))

    for allele in "pqr":
        print("{}:\t{}\t{}%".format(allele, alleleDict[allele], 100*alleleDict[allele]/sum(alleleDict.values())))
    '''
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
    freqTable = [["A", "I", "O"]]
    for gen in range(nGen):
        pop = newGen(pop)
        phenFreq = calcPhenoFreq(pop)
        for ind in pop[0]:
            ind.calcFec(phenFreq)
        for ind in pop[1]:
            ind.calcFec(phenFreq)
        freqTable.append([str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"])])
    genoFreq = calcGenoFreq(pop)
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'w') as outfile:
            for line in freqTable:
                print("\t".join(line), file=outfile)
    else:
        for line in freqTable:
            print("\t".join(line))
    #print(genoFreq)
