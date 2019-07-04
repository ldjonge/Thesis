#!/usr/bin/env python3
import random
from standardFunctions import *
from numpy.random import choice
import numpy as np
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
    data = pd.read_csv(file, sep=",", header=None)
    data = data.div(data.sum(axis=1), axis=0)
    return(data)

# Define Males
class Male:
    # New individuals are created as offspring from two parental alleles, in males genotype has no effect on phenotype
    def __init__(self, pAll, mAll):
        self.genotype = pAll+mAll
        self.phenotype = "M"
        self.taken = 0
        self.migrate = 0

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
                    self.fecundity = popInfo["Ofec"]
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
                    self.fecundity = popInfo["Ifec"]
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
                    self.fecundity = popInfo["Afec"]
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

    def eggLay(self, pop, Eggs):
        nEggs = Eggs * self.fecundity
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

def startingPop(popInfo, params):
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
        ind.learning(totalPop, params)
    for ind in femalePop:
        ind.calcFec(popInfo)
    pop = [malePop, femalePop]
    return pop

# Phenotype distribution is relevant in many cases, sometimes male frequency is also important but in most cases it is not
def calcPhenoFreq(pop, male=False):
    freqDict = {"M":0, "A":0, "I":0, "O":0}
    if type(pop[0]) == list:
        for ind in pop[0]:
            freqDict[ind.phenotype] += 1
        for ind in pop[1]:
            freqDict[ind.phenotype] += 1
    else:
        for ind in pop:
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
def matingSearch(pop, params, popDict):
    # The number of interactions per population/generation should be recorded
    matings = 0
    contacts = 0
    MMcontacts = 0
    deaths = 0
    totalPop = pop[0] + pop[1]
    N = len(totalPop)
    if len(pop[1])!=0:
        for male in pop[0]:
            if male.taken == 0:
                # the chance of approaching a specific individual is determined for each potential pair, based on the male's preference.
                hitChance = []
                for ind in pop[0]:
                    # There should be the potential for a male to approach another male
                    if ind.taken != 0:
                        hitChance.append(0)
                    else:
                        hitChance.append(max(params["maleRec"]*male.aPref, 0))
                for fem in pop[1]:
                    if fem.taken != 0:
                        hitChance.append(0)
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
                    if random.random() < male.mSucc*N/popDict["K"]:    # Male mating success is used as probability of finding a mate
                        mate = choice(totalPop, p=hitChance)
                    else:
                        mate = None
                    if type(mate)==Female:
                        contacts += 1
                        # Female mating success incorporates both the chance of contact leading to mating and the chance of mating leading to fertilisation
                        if random.random() <= mate.mSucc*male.fertility:
                            matings += 1
                            mate.mate(male, params)
                            mate.mSucc *= params["mateFFertEff"]
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
        deaths += (N - len(pop[0]) - len(pop[1]))
    else:
        pass
    return [matings, contacts, MMcontacts, deaths]

# The new population size will be decided based on the number of eggs produced in the previous generation and  the carrying capacity of the population
def newPopSize(nEggs, pop, K):
    oldPopSize = len(pop)
    avgPop = 4*oldPopSize/nEggs
    newPopSize = randomRound(avgPop*(np.exp(0.5*(K-avgPop)/K)))
    if avgPop > K:
        newPopSize = max(K, newPopSize)
    return(newPopSize)

# individuals will be chosen randomly from the population of eggs
def popControl(pop, size):
    newPop = list(choice(pop, size=int(size), replace=False))
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
# These parameters are recorded per sex, as well as per female morph
def recordFecStats(pop):
    if len(pop[0]) > 0:
        totalMFer = 0
        totalMMSucc = 0
        totalMSurv = 0
        for m in pop[0]:
            totalMFer += m.fertility
            totalMMSucc += m.mSucc
            totalMSurv += m.surv
        avgMFer = totalMFer/len(pop[0])
        avgMMSucc = totalMMSucc/len(pop[0])
        avgMSurv = totalMSurv/len(pop[0])
    else:
        avgMFer = np.nan
        avgMMSucc = np.nan
        avgMSurv = np.nan
    if len(pop[1]) > 0:
        totalFFec = 0
        totalFMSucc = 0
        totalFSurv = 0
        totalAFec = 0
        totalIFec = 0
        totalOFec = 0
        totalAMSucc = 0
        totalIMSucc = 0
        totalOMSucc =0
        totalASurv = 0
        totalISurv = 0
        totalOSurv = 0
        totalA = 0
        totalI = 0
        totalO = 0
        for f in pop[1]:
            totalFFec += f.fecundity
            totalFMSucc += f.mSucc
            totalFSurv += f.surv
            if f.phenotype == "A":
                totalAFec += f.fecundity
                totalAMSucc += f.mSucc
                totalASurv += f.surv
                totalA += 1
            elif f.phenotype == "I":
                totalIFec += f.fecundity
                totalIMSucc += f.mSucc
                totalISurv += f.surv
                totalI += 1
            elif f.phenotype == "O":
                totalOFec += f.fecundity
                totalOMSucc += f.mSucc
                totalOSurv += f.surv
                totalO += 1
        avgFFec = totalFFec/len(pop[1])
        avgFMSucc = totalFMSucc/len(pop[1])
        avgFSurv = totalFSurv/len(pop[1])
        if totalA > 0:
            avgAFec = totalAFec/totalA
            avgAMSucc = totalAMSucc/totalA
            avgASurv = totalASurv/totalA
        else:
            avgAFec = np.nan
            avgAMSucc = np.nan
            avgASurv = np.nan
        if totalI > 0:
            avgIFec = totalIFec/totalI
            avgIMSucc = totalIMSucc/totalI
            avgISurv = totalISurv/totalI
        else:
            avgIFec = np.nan
            avgIMSucc = np.nan
            avgISurv = np.nan
        if totalO > 0:
            avgOFec = totalOFec/totalO
            avgOMSucc = totalOMSucc/totalO
            avgOSurv = totalOSurv/totalO
        else:
            avgOFec = np.nan
            avgOMSucc = np.nan
            avgOSurv = np.nan
    else:
        avgFFec = np.nan
        avgFMSucc = np.nan
        avgFSurv = np.nan
        avgAFec = np.nan
        avgAMSucc = np.nan
        avgASurv = np.nan
        avgIFec = np.nan
        avgIMSucc = np.nan
        avgISurv = np.nan
        avgOFec = np.nan
        avgOMSucc = np.nan
        avgOSurv = np.nan
    return [avgMFer, avgMMSucc, avgMSurv, avgFFec, avgFMSucc, avgFSurv, avgAFec, avgAMSucc, avgASurv, avgIFec, avgIMSucc, avgISurv, avgOFec, avgOMSucc, avgOSurv]

# Preference should be recorded, again as a population average
def recordPref(pop):
    totalAPref = 0
    totalIPref =0
    totalOPref = 0
    for m in pop[0]:
        totalAPref += m.aPref
        totalIPref += m.iPref
        totalOPref += m.oPref
    prefSum = totalAPref + totalIPref + totalOPref
    if prefSum > 0:
        avgAPref = totalAPref/prefSum
        avgIPref = totalIPref/prefSum
        avgOPref = totalOPref/prefSum
        return [avgAPref, avgIPref, avgOPref]
    else:
        return [np.nan,np.nan,np.nan]

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
				migOut = random.random()
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
	return migrationList

def runSim():
    paramDict = readParams()
    params = readPopInfo()
    print(len(params))
    migrationMatrix = readMigration()
    pops = []
    nPop = min(len(params), len(migrationMatrix.index))
    print(nPop)
    for i in range(nPop):
        pop = params[i]
        pops.append(startingPop(pop, paramDict))
    nEggs = paramDict["nEggs"]
    freqTable = [["Pop", "Gen", "Pheno", "Value"]]
    for gen in range(1, int(paramDict["nGen"])+1):#First generation should be called 1
        totalLen = 0
        matings = [0 for i in range(nPop)]
        contacts = [0 for i in range(nPop)]
        MMcontacts = [0 for i in range(nPop)]
        migrations = [0 for i in range(nPop)]
        deaths = [0 for i in range(nPop)]
        popSizes = []
        for pop in pops:
            id = pops.index(pop)+1
            prefs = recordPref(pop)
            phenFreq = calcPhenoFreq(pop)
            freqTable.append([id, gen, "preAPref", prefs[0]])
            freqTable.append([id, gen, "preIPref", prefs[1]])
            freqTable.append([id, gen, "preOPref", prefs[2]])
            freqTable.append([id, gen, "A", phenFreq["A"]])
            freqTable.append([id, gen, "I", phenFreq["I"]])
            freqTable.append([id, gen, "O", phenFreq["O"]])
            freqTable.append([id, gen, "M", len(pop[0])])
            freqTable.append([id, gen, "F", len(pop[1])])
            freqTable.append([id, gen, "T", len(pop[0])+len(pop[1])])
            popSizes.append(len(pop[0])+len(pop[1]))
        for i in range(int(paramDict["genLength"])):
            for pop in pops:
                id = pops.index(pop)
                results = matingSearch(pop, paramDict, params[id])
                matings[id] += results[0]
                contacts[id] += results[1]
                MMcontacts[id] += results[2]
                deaths[id] += results[3]
                for fem in pop[1]:
                    if fem.taken != 0:
                        fem.taken -= 1
            cyclMigr = migrate(pops, migrationMatrix)
            migrations = [a+b for a,b in zip(migrations, cyclMigr)]
        newPops = []
        totalLen = 0
        for pop in pops:
            id = pops.index(pop)
            avgFecs = recordFecStats(pop)
            prefs = recordPref(pop)
            freqTable.append([id+1, gen, "MFer", avgFecs[0]])
            freqTable.append([id+1, gen, "MMSucc", avgFecs[1]])
            freqTable.append([id+1, gen, "MSurv", avgFecs[2]])
            freqTable.append([id+1, gen, "FFec", avgFecs[3]])
            freqTable.append([id+1, gen, 'FMSucc', avgFecs[4]])
            freqTable.append([id+1, gen, 'FSurv', avgFecs[5]])
            freqTable.append([id+1, gen, 'AFec', avgFecs[6]])
            freqTable.append([id+1, gen, 'AMSucc', avgFecs[7]])
            freqTable.append([id+1, gen, 'ASurv', avgFecs[8]])
            freqTable.append([id+1, gen, 'IFec', avgFecs[9]])
            freqTable.append([id+1, gen, 'IMSucc', avgFecs[10]])
            freqTable.append([id+1, gen, 'ISurv', avgFecs[11]])
            freqTable.append([id+1, gen, 'OFec', avgFecs[12]])
            freqTable.append([id+1, gen, 'OMSucc', avgFecs[13]])
            freqTable.append([id+1, gen, 'OSurv', avgFecs[14]])
            freqTable.append([id+1, gen, 'Matings', matings[id]])
            freqTable.append([id+1, gen, 'Contacts', contacts[id]])
            freqTable.append([id+1, gen, 'MMContacts', MMcontacts[id]])
            freqTable.append([id+1, gen, 'Migrations', migrations[id]])
            freqTable.append([id+1, gen, 'Deaths', deaths[id]])
            freqTable.append([id+1, gen, "APref", prefs[0]])
            freqTable.append([id+1, gen, "IPref", prefs[1]])
            freqTable.append([id+1, gen, "OPref", prefs[2]])
            newPop = []
            for f in pop[1]:
                f.eggLay(newPop, paramDict["nEggs"])
            size = newPopSize(paramDict["nEggs"], newPop, params[id]["K"])
            newPop = popControl(newPop, size)
            popSize = len(newPop[0]) + len(newPop[1])
            if popSize > 0 and popSizes[id] ==0:
                print("New Population formed at {} in generation {}".format(id+1, gen))
            elif popSize == 0 and popSizes[id] > 0:
                print("Population {} extinct in generation {}".format(id+1, gen))
            totalLen += popSize
            for ind in newPop[0]:
                ind.calcFec(popInfo=params[id])
                ind.learning(newPop[0]+newPop[1], paramDict)
            for ind in newPop[1]:
                ind.calcFec(popInfo=params[id])
            newPops.append(newPop)
        if gen%10 == 0:
            print("Generation {} complete, population size {}".format(str(gen), str(totalLen)))
        pops = newPops
    return freqTable

if __name__ == "__main__":
    freqTable = pd.DataFrame(runSim())
    if len(sys.argv) > 1:
        freqTable.to_csv(sys.argv[1], header=None, index=None)
