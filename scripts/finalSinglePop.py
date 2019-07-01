#!/usr/bin/env python3
import random
from standardFunctions import *
from numpy.random import choice
from numpy import exp
import sys
import math
import pandas as pd
from copy import copy

# Create a dictionary for all possible genotypes to refer to a phenotype
phenoDict = {"pp":"A", "pq":"A", "pr":"A", "qp":"A", "rp":"A", "qq":"I", "qr":"I", "rq":"I", "rr":"O"}

# Function to read in parameter values from external file, allowing for a potential alternative file
# The function returns a dictionary with the values for each parameter
def readParams(file="paramFiles/parameters.ini"):
    with open(file, "r") as infile:
        paramDict = {}
        for line in infile:
            line = line.strip()
            line = line.split(" ") # ini file should be space separated
            # Numeric parameters should be converted to floats, other parameters should stay strings.
            # In case integers are needed for a specific parameter this conversion can be performed later
            try:
                # Numeric parameters should be converted to floats, other parameters should stay strings.
                # In case integers are needed for a specific parameter this conversion can be performed later
                paramDict[line[0]] = float(line[1])
            except ValueError:
                paramDict[line[0]] = line[1]
    return(paramDict)

def readPopInfo(file="paramFiles/popInfo.csv"):
    with open(file, 'r') as data:
        header = data.readline()
        header = header.strip()
        header = header.split(",")
        pops = []
        for line in data.readlines():
            line = line.strip()
            line = line.split(",")
            popDict={}
            for i in range(len(header)):
                popDict[header[i]] = float(line[i])
            pops.append(popDict)
        return pops

def readMigration(file="paramFiles/dispersalMatrix.csv"):
    data = pandas.read_csv(file, sep=",", header=None)
    data = data.div(data.sum(axis=1), axis=0)
    return(data)

# Define Males
class Male:
    # New individuals are created as offspring from two parental alleles, in males genotype has no effect on phenotype
    def __init__(self, pAll, mAll):
        self.genotype = pAll+mAll
        self.phenotype = "M"
        self.taken = 0

    # Mating Success and Survival chance are specified in an external file, however if these are not present baseline values of 1 will be applied.
    def calcFec(self, popInfo=None):
        if popInfo == None:
            self.fertility = 1
            self.mSucc = 1
            self.surv = 1
        else:
            try:
                self.fertility = popInfo["MFert"]
            except KeyError:
                self.fertility = 1
            try:
                self.mSucc= popInfo["MMSucc"]
            except KeyError:
                self.mSucc = 1
            try:
                self.surv = popInfo["MSurv"]
            except KeyError:
                self.surv = 1

    # Males go through a learning process based on the frequencies of morphs in the population around them
    def learning(self, pop, params):
        if len(pop) >= 50:
            selected = choice(pop, size=50, replace=False)
        else:
            selected = pop
        prefDict = calcPhenoFreq(selected, male=True)
        self.aPref = params["mPref"]*prefDict["M"]+params["aPref"]*prefDict["A"]
        self.iPref = params["iPref"]*prefDict["I"]+params["hPref"]
        self.oPref = params["oPref"]*prefDict["O"]+params["hPref"]

    # In case a full population should be shown males will be represented by "M"
    def __str__(self):
        return self.phenotype

# Define Females
class Female:
    # New individuals are created from two parental alleles, genotype determines female phenotype.
    def __init__(self,pAll,mAll):
        self.genotype = pAll+mAll
        self.phenotype = phenoDict[self.genotype]
        self.taken = 0 # A mating female is taken out of the mating search for a short time
        self.mates = [] # A list to store the individuals a female has mated with for genotype determination of offspring
        self.migrate = 0 # A counter to keep track of migration

    # Fecundity, Mating Success and Survival chance are specified in an external file, however if these are not present baseline values of 1 will be applied.
    def calcFec(self, popInfo=None):
        if popInfo == None:
            self.fecundity = 1
            self.mSucc = 1
            self.surv = 1
        else:
            if self.phenotype == "O":
                try:
                    self.fertility = popInfo["Ofec"]
                except KeyError:
                    self.fecundity = 1
                try:
                    self.mSucc = popInfo["OMSucc"]
                except KeyError:
                    self.mSucc = 1
                try:
                    self.surv = popInfo["OSurv"]
                except KeyError:
                    self.surv = 1
            elif self.phenotype == "I":
                try:
                    self.fertility = popInfo["Ifec"]
                except KeyError:
                    self.fecundity = 1
                try:
                    self.mSucc = popInfo["IMSucc"]
                except KeyError:
                    self.mSucc = 1
                try:
                    self.surv = popInfo["ISurv"]
                except KeyError:
                    self.surv = 1
            elif self.phenotype == "A":
                try:
                    self.fertility = popInfo["Afec"]
                except KeyError:
                    self.fecundity = 1
                try:
                    self.mSucc = popInfo["AMSucc"]
                except KeyError:
                    self.mSucc = 1
                try:
                    self.surv = popInfo["ASurv"]
                except KeyError:
                    self.surv = 1

    # When a female mates, fecundity is affected, and she is removed from mating searches for two cycles. The male is added to the list of males to determine the distribution of the offspring
    def mate(self, male, params):
        self.fecundity *= params["mateFecEff"]
        self.mates.append(male.genotype)
        self.taken += 2 #Copulation is assumed to take up 2 additional cycles worth of time
        male.taken += 2

    # Fertilised females are assumed to lay all eggs at the same time, with the sperm used for fertilisation depending on the order of mating
    # Genotypes of the offspring are evenly distributed according to Mendelian genetics, with randomness of survival accounting for the realistic randomness of genotype distribution in adults

    def eggLay(self, pop, nEggs):
        nEggs = nEggs * self.fecundity
        nMates = len(self.mates)
        if nMates > 0:
            for male in self.mates[:-1]:
                for i in range(randomRound(0.265/(nMates-1)*nEggs/8)):
                    pop.append(Male(male[0], self.genotype[0]))
                    pop.append(Male(male[0], self.genotype[1]))
                    pop.append(Male(male[1], self.genotype[0]))
                    pop.append(Male(male[1], self.genotype[1]))
                    pop.append(Female(male[0], self.genotype[0]))
                    pop.append(Female(male[0], self.genotype[1]))
                    pop.append(Female(male[1], self.genotype[0]))
                    pop.append(Female(male[1], self.genotype[1]))
            male = self.mates[-1]
            for i in range(randomRound(0.735*nEggs/8)):
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
            for i in range(randomRound(nEggs/8)):
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

    # In case a full population should be shown females will be represented by their phenotype
    def __str__(self):
        return self.phenotype

# Genotypes of a starting population are determined randomly based on allele frequencies
def randomAllele(p,q,r):
    newRand = random.random()
    if newRand <= p:
        allele = "p"
    elif newRand <= p+q:
        allele = "q"
    else:
        allele = "r"
    return allele

# Male and female starting populations are created separately
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

# Phenotype distribution is relevant in many cases, sometimes male frequency is also important but in most cases it is not
def calcPhenoFreq(pop, male=False):
    freqDict = {"M":0, "A":0, "I":0, "O":0}
    for ind in pop[0]:
        freqDict[ind.phenotype] += 1
    for ind in pop[1]:
        freqDict[ind.phenotype] += 1
    phenoDist = {}
    if male==True:
        for phen in "MAIO":
            try:
                phenoDist[phen] = freqDict[phen]/sum([freqDict["A"], freqDict["I"], freqDict["O"], freqDict["M"]])
            except ZeroDivisionError:
                phenoDist[phen] = 0
    else:
        for phen in "AIO":
            try:
                phenoDist[phen] = freqDict[phen]/sum([freqDict["A"], freqDict["I"], freqDict["O"]])
            except ZeroDivisionError:
                phenoDist[phen] = 0
    return phenoDist

# Males will search for mates based on their preference.
def matingSearch(pop, K, params):
    # The number of interactions per population/generation should be recorded
    matings = 0
    contacts = 0
    MMcontacts = 0
    if len(pop[1])!=0:
        for male in pop[0]:
            if male.taken == 0:
                # the chance of approaching a specific individual is determined for each potential pair, based on the male's preference.
                hitChance = []
                for ind in pop[0]:
                    # There should be the potential for a male to approach another male
                    hitChance.append(max(params["maleRec"]*male.aPref, 0))
                for fem in pop[1]:
                    if fem.taken != 0:
                        pass
                    elif fem.phenotype == "A":
                        hitChance.append(max(male.aPref, 0))
                    elif fem.phenotype == "I":
                        hitChance.append(max(male.iPref, 0))
                    elif fem.phenotype == "O":
                        hitChance.append(max(male.oPref, 0))
                totalW = sum(hitChance)
                # Probabilities of selecting a specific individual are calculated from the weights
                if totalW != 0:
                    for i in range(len(hitChance)):
                        hitChance[i] = hitChance[i]/totalW
                    if random.random() < male.mSucc*N/K:    # Male mating success is used as probability of finding a mate
                        mate = choice(totalPop, p=hitChance)
                    else:
                        mate = None
                    if type(mate)==Female:
                        contacts += 1
                        # Female mating success incorporates both the chance of contact leading to mating and the chance of mating leading to fertilisation
                        if random.random() <= mate.mSucc*male.fertility:
                            matings += 1
                            mate.mate(male)
                            mate.mSucc *= params["mateFertEff"]
                            male.fertility *= params["mateMFertEff"]
                            mate.surv *= params["mateSurvEff"]
                            male.surv *= params["mateSurvEff"]
                            if mate.phenotype == "A":
                                male.aPref *=params["matePrefEff"]
                            elif mate.phenotype == "I":
                                male.iPref *=params["matePrefEff"]
                            elif mate.phenotype == "O":
                                male.oPref *=params["matePrefEff"]
                        else:
                            mate.mSucc *=params["failFertEff"]
                            mate.surv *= params["failSurvEff"]


                            if mate.phenotype == "A":
                                male.aPref *=params["failPrefEff"]
                            elif mate.phenotype == "I":
                                male.iPref *=params["failPrefEff"]
                            elif mate.phenotype == "O":
                                male.oPref *=params["failPrefEff"]
                    elif type(mate) == Male:
                        MMcontacts += 1
                        mate.fertility *= params["failFertEff"]
                        mate.surv *= params["failSurvEff"]
                        male.aPref *=params["failPrefEff"]
            else:
                male.taken -= 1
        # Individuals that die will be removed from the population
        pop[0] = [i for i in pop[0] if i.surv > random.random()]
        pop[1] = [i for i in pop[1] if i.surv > random.random()]
    else:
        pass
    return [matings, contacts, MMcontacts]

# The new population size will be decided based on the number of eggs produced in the previous generation and  the carrying capacity of the population
def newPopSize(nEggs, pop, K):
    oldPopSize = len(pop)
    avgPop = 2*oldPopSize/nEggs
    newPopSize = randomRound(avgPop*(exp(0.5*(K-avgPop)/K)))
    if avgPop > K:
        newPopSize = max(K, newPopSize)
    return(newPopSize)

# individuals will be chosen randomly from the population of eggs
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
    newPop = [malePop, femalePop]
    return(newPop)

# Fecundity parameters should be recorded, not per individual but as a population average


def migrate(pops, migrationMatrix):
	newPops = [([],[]) for i in pops]
	migrationList = []
	for pop in pops:
		migrations = 0
		id = pops.index(pop)
		mRates = migrationMatrix.iloc[id]
		mChance = mRates.cumsum()
		for male in pop[0]:
			if male.migrate == 0:	#Individuals that are still migrating cannot migrate again
				migChance = random.random()
				for i in range(len(mChance)):
					if migOut <= mChance[i]:
						newPop = i
						break
				if newPop != id:
					migrations += 1
					male.migrate += random.randint(1,5)	#Depending on distance and speed migration may take a variable amount of time
				newPops[newPop][0].append(male)	#Individuals are added to the population for the next cycle
			else:
				male.migrate -= 1	#Migrating individuals get one step closer to their destination
				newPops[id][0].append(male)	#But are staying in the population they are going towards
		for female in pop[1]:
			if female.migrate == 0:
				migChance = random.random()
				for i in range(len(mChance)):
					if migOut <= mChance[i]:
						newPop = i
						break
				if newPop != id:
					migrations += 1
					female.migrate += 2
				newPops[newPop][0].append(female)
			else:
				female.migrate -= 1
				newPops[id][0].append(female)
		migrationList.append(migrations)
	pops = newPops
	return migrations

def runSim():
	paramDict = readParams()
	params = readPopInfo()
	migrationMatrix = readMigration()
	pops = []
	nPop = min(len(params), len(migrationMatrix.index))
	for i in range(nPop):
		pop = params[i]
		pops.append(startingPop(pop))
	nEggs = paramDict["nEggs"]
	freqTable = [["Pop", "Gen", "Pheno", "Value"]]
	for gen in range(1, int(paramDict["nGen"])+1):	#First generation should be called 1
		totalLen = 0
		matings = [0 for i in range(nPop)]
		contacts = [0 for i in range(nPop)]
		MMcontacts = [0 for i in range(nPop)]
		migrations = [0 for i in range(nPop)]
		for pop in pops:
			id = pops.index(pop)+1
			prefs = recordPref(pop)
			freqTable.append([id, gen, "preAPref", prefs[0]])
			freqTable.append([id, gen, "preIPref", prefs[1]])
			freqTable.append([id, gen, "preOPref", prefs[2]])
		for i in range(paramDict["genLength"]):
			for pop in pops:
				id = pops.index(pop)
				results = matingSearch(pop, """Whatever else is needed here""")
				matings[id] += results[0]
				contacts[id] += results[1]
				MMcontacts[id] += results[2]
				for fem in pop[1]:
				if fem.taken != 0:
					fem.taken -= 1
			cyclMigr = migrate(pops, migrationMatrix)
			migrations = [a+b for a,b in zip(migrations, cyclMigr)]
		for pop in pops:
			id = pops.index(pop)
			avgFecs = recordFec(pop)
			prefs = recordPref(pop)
			freqTable.append("""STUFF, lets figure out what the record functions spit out first""")
