#!/usr/bin/env python3
import random
from standardFunctions import *
from numpy.random import choice
import sys

phenoDict = {"pp":"A", "pq":"A", "pr":"A", "qp":"A", "rp":"A", "qq":"I", "qr":"I", "rq":"I", "rr":"O"}

#Starting parameters
N = 200
p=2/6
q=2/6
r=1/3
nGen = 100

class Male:
    def __init__(self, pAll, mAll):
        self.genotype = pAll+mAll
        self.phenotype = "M"
        self.fecundity = 1

class Female:
    def __init__(self,pAll,mAll, phenoFreq):
        self.genotype = pAll+mAll
        self.phenotype = phenoDict[self.genotype]
        self.fecundity = calcNFDF(self.phenotype, phenoFreq)

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

def createFemalePop(N, p, q, r, phenoFreq):
    population = []
    pAll = randomAllele(p,q,r)
    mAll = randomAllele(p,q,r)
    for i in range(N):
        pAll = randomAllele(p,q,r)
        mAll = randomAllele(p,q,r)
        population.append(Female(pAll,mAll, phenoFreq))
    return(population)

def reproduce(male, female, phenoFreq):
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
        return Female(pAll, mAll, phenoFreq)

def newGen(pop):
    newMalePop = []
    newFemalePop = []
    phenoFreq = calcPhenoFreq(pop)
    maleFec = [ind.fecundity for ind in pop[0]]
    femaleFec = [ind.fecundity for ind in pop[1]]
    for i in range(N):
        pat = choice(pop[0], 1, maleFec)[0]
        mat = choice(pop[1], 1, femaleFec)[0]
        offspring = reproduce(pat, mat, phenoFreq)
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
    phenFreq = phenoFreq(p,q,r)
    malePop = createMalePop(N//2, p,q,r)
    femalePop = createFemalePop(N//2, p,q,r, phenFreq)
    pop = (malePop, femalePop)
    freqTable = [["A", "I", "O"]]
    for gen in range(nGen):
        pop = newGen(pop)
        phenFreq = calcPhenoFreq(pop)
        freqTable.append([str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"])])
    genoFreq = calcGenoFreq(pop)
    with open(sys.argv[1], 'w') as outfile:
        for line in freqTable:
            print("\t".join(line), file=outfile)
    #print(genoFreq)
