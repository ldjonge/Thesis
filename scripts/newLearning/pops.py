import numpy as np
from baseFunctions import *


# Create a dictionary for all possible genotypes to refer to a phenotype
phenoDict = {"pp":"A", "pq":"A", "pr":"A", "qp":"A", "rp":"A", "qq":"I", "qr":"I", "rq":"I", "rr":"O"}

class Male:
    def __init__(self, pAll, mAll, params):
        self.genotype = pAll+mAll # New individuals are created as offspring from two parental alleles
        self.phenotype = "M"    # In males genotype has no effect on phenotype
        self.taken = 0
        self.migrate = 0
        self.prefs = {"A": params["{}prefA".format(self.phenotype)], "I": params["{}prefI".format(self.phenotype)], "O": params["{}prefO".format(self.phenotype)]}

    # Fertility and Survival chance are specified in an external file, however if these are not present baseline values of 1 will be applied.
    def calcFec(self, popInfo=None):
        if popInfo == None:
            self.fertility = 1
            self.surv = 1
        else:
            try:
                self.fertility = popInfo["MFert"]
            except KeyError:
                self.fertility = 1
            try:
                self.surv = popInfo["MSurv"]
            except KeyError:
                self.surv = 1

    def calcmSucc(self):
        valsum = sum(self.prefs.values())
        #for key in self.prefs.keys():
        #    self.prefs[key] = (self.prefs[key]/valsum)*0.7+0.1
        #if sum(self.prefs.values()) < 0.99 or sum(self.prefs.values()) > 1.01:
        #    print(self.prefs)
        self.mSucc = (max(self.prefs.values())/(sum(self.prefs.values())))

    # Males go through a learning process based on the frequencies of morphs in the population around them
    def learning(self, pop, params):
        # If learning pre-mating:
        chance = [fem.vis for fem in pop]
        chSum = sum(chance)
        chance = [i/chSum for i in chance]
        for fem in np.random.choice(pop, size=min(len(pop), 10), p=chance, replace=False):
            # if random.random() < fem.mSucc:
            self.prefs[fem.phenotype] *= params["matePrefEff"]**0.5
            #else:
            #    self.prefs[fem.phenotype] *= params["failPrefEff"]

    # In case a full population should be shown males will be represented by "M"
    def __str__(self):
        return self.phenotype

class Female:
    def __init__(self,pAll,mAll):
        self.genotype = pAll+mAll   # New individuals are created from two parental alleles
        self.phenotype = phenoDict[self.genotype]   # Genotype determines female phenotype.
        self.taken = 0 # A mating female is taken out of the mating search for a short time
        self.mates = [] # A list to store the individuals a female has mated with for genotype determination of offspring
        self.migrate = 0 # A counter to keep track of migration (later?)

    # Fecundity, Mating Success and Survival chance are specified in an external file, however if these are not present baseline values of 1 will be applied.
    def calcFec(self, popInfo=None):
        if popInfo == None:
            self.fecundity = 400
            self.mSucc = 1
            self.surv = 1
        else:
            try:
                self.fecundity = popInfo["{}fec".format(self.phenotype)]
            except KeyError:
                self.fecundity = 400
            try:
                self.mSucc = popInfo["{}MSucc".format(self.phenotype)]
            except KeyError:
                self.mSucc = 1
            try:
                self.surv = popInfo["{}Surv".format(self.phenotype)]
            except KeyError:
                self.surv = 1

    def calcVis(self, phenFreq):
        if self.phenotype == "A":
            self.vis = 1 #phenFreq["A"]**0.5
        else:
            self.vis = 1
    """
    When a female mates, fecundity is affected, and she is removed from mating searches for two cycles.
    The male is added to the list of males to determine the distribution of the offspring, and also removed from mating searches for two cycles.
    """
    def mate(self, male, params):
        self.fecundity *= params["mateFecEff"]
        if random.random() <= male.fertility: # Even if a male does not successfully fertilise a female, the rest of the mating process still occurs.
            self.mates.append(male.genotype)
        self.taken += 2 #Copulation is assumed to take up 2 additional cycles worth of time
        male.taken += 2

    """Fertilised females are assumed to lay all eggs at the same time, with the sperm used for fertilisation
    depending on the order of mating. Genotypes of the offspring are evenly distributed according to Mendelian
    genetics, with randomness of survival accounting for the realistic randomness of genotype distribution in adults
    """
    def eggLay(self, pop, params):
        nEggs = np.random.normal(1, 0.1)*self.fecundity
        nMates = len(self.mates)
        if nMates > 1:
            for male in self.mates[:-1]:
                for i in range(randomRound(0.265/(nMates-1)*nEggs/8)):
                    pop.append(Male(male[0], self.genotype[0], params))
                    pop.append(Male(male[0], self.genotype[1], params))
                    pop.append(Male(male[1], self.genotype[0], params))
                    pop.append(Male(male[1], self.genotype[1], params))
                    pop.append(Female(male[0], self.genotype[0]))
                    pop.append(Female(male[0], self.genotype[1]))
                    pop.append(Female(male[1], self.genotype[0]))
                    pop.append(Female(male[1], self.genotype[1]))
            male = self.mates[-1]
            for i in range(randomRound(0.735*nEggs/8)):
                pop.append(Male(male[0], self.genotype[0], params))
                pop.append(Male(male[0], self.genotype[1], params))
                pop.append(Male(male[1], self.genotype[0], params))
                pop.append(Male(male[1], self.genotype[1], params))
                pop.append(Female(male[0], self.genotype[0]))
                pop.append(Female(male[0], self.genotype[1]))
                pop.append(Female(male[1], self.genotype[0]))
                pop.append(Female(male[1], self.genotype[1]))
        elif nMates == 1:
            male = self.mates[-1]
            for i in range(randomRound(nEggs/8)):
                pop.append(Male(male[0], self.genotype[0], params))
                pop.append(Male(male[0], self.genotype[1], params))
                pop.append(Male(male[1], self.genotype[0], params))
                pop.append(Male(male[1], self.genotype[1], params))
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

"""
Phenotype distribution is relevant in many cases.
While male frequency may be relevant in some cases, it usually is not, hence the default value
"""
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
        for phen in "MA":
            try:
                phenoDist[phen] = freqDict[phen]/sum([freqDict["A"], freqDict["M"]])
            except ZeroDivisionError:
                phenoDist[phen] = 0
    else:
        for phen in "AIO":
            try:
                phenoDist[phen] = freqDict[phen]/sum([freqDict["A"], freqDict["I"], freqDict["O"]])
            except ZeroDivisionError:
                phenoDist[phen] = 0
    return phenoDist

# Male and female starting populations are created separately
def createMalePop(N, p, q, r, params):
    population = []
    for i in range(N):
        pAll = randomAllele(p,q,r)
        mAll = randomAllele(p,q,r)
        population.append(Male(pAll,mAll, params))
    return(population)

def createFemalePop(N, p, q, r):
    population = []
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
    malePop = createMalePop(N//2, p,q,r, popInfo)
    femalePop = createFemalePop(N//2, p,q,r)
    phenFreq = calcPhenoFreq((malePop, femalePop), male=True)
    totalPop = malePop+femalePop
    for ind in femalePop:
        ind.calcFec(popInfo)
        ind.calcVis(phenFreq)
    for ind in malePop:
        ind.calcFec(popInfo)
        #ind.learning(femalePop, params)
        ind.calcmSucc()
    pop = [malePop, femalePop]
    return pop

def newPopSize(pop, K):
    oldPopSize = len(pop)
    avgPop = oldPopSize/100 # 1% of eggs are expected to survive on average
    newPopSize = randomRound(avgPop*(np.exp(0.5*(K-avgPop)/K))) #Adjust the number of survivors based on carrying capacity
    """
    At large numbers of eggs the formula will reduce population size below carrying capacity
    which is obviously unrealistic. While this seems like something that needs a change it may not be needed
    due to the low chance of these numbers actually occurring
    """
    if avgPop > K:
        newPopSize = max(K, newPopSize)
    """
    The number off survivors should be randomised, as situations will vary over the years
    """
    #newPopSize = randomRound(np.random.normal(1, 0.25)*newPopSize)
    return(newPopSize)

# Individuals will be chosen randomly from the population of eggs
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
