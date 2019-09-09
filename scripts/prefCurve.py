import random
import numpy as np
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

class Male:
    def __init__(self, pAll, mAll):
        self.genotype = pAll+mAll # New individuals are created as offspring from two parental alleles
        self.phenotype = "M"    # In males genotype has no effect on phenotype
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
        self.prefs = {"A": params["{}prefA".format(self.phenotype)], "I": params["{}prefI".format(self.phenotype)], "O": params["{}prefO".format(self.phenotype)]}
        # If learning pre-mating:
        for fem in np.random.choice(pop, size=min(len(pop), 10), replace=False):
            if random.random() < fem.mSucc:
                self.prefs[fem.phenotype] *= params["succEff"]
            else:
                self.prefs[fem.phenotype] *= params["failEff"]

    # In case a full population should be shown males will be represented by "M"
    def __str__(self):
        return self.phenotype

class testFem:
    def __init__(self, phen):
        self.phenotype = phen
        if phen == "A":
            self.mSucc = 0.8
        else:
            self.mSucc = 0.8

    def __str__(self):
        return self.phenotype

params = {"MprefA": 1, "MprefI": 1, "MprefO": 1, "succEff": 3, "failEff": 0.9}

A = range(51)
options = [v for v in itertools.product(A, repeat=3) if sum(v)==50]



def testPref(malePop, testPop, params):
    totalPref = {"A": 0, "I": 0, "O": 0}
    for male in malePop:
        male.learning(testPop, params)
        avg = sum(male.prefs.values())
        for i in totalPref.keys():
            totalPref[i] += (male.prefs[i]/avg)
    return totalPref

outTable = {"A": [], "I": [], "O": [], "Apref": [], "Ipref": [], "Opref": []}
for i in options:
    freqs = i[0]*"A"+i[1]*"I"+i[2]*"O"
    outTable["A"].append(i[0])
    outTable["I"].append(i[1])
    outTable["O"].append(i[2])
    testPop = [testFem(phen) for phen in freqs]
    malePop = [Male("p", "p") for i in range(100)]
    pref = testPref(malePop, testPop, params)
    outTable["Apref"].append(pref["A"])
    outTable["Ipref"].append(pref["I"])
    outTable["Opref"].append(pref["O"])

out = pd.DataFrame.from_dict(outTable)
fig, axarr = plt.subplots(1,1, sharex =False, sharey=False, figsize=(10,8))
plot1 = sns.scatterplot("A", "Apref", data=out, ax=axarr)
fig.savefig("nprefTestA.png")
fig, axarr = plt.subplots(1,1, sharex =False, sharey=False, figsize=(10,8))
plot3 = sns.scatterplot("O", "Opref", data=out, ax=axarr)
fig.savefig("nprefTestO.png")
fig, axarr = plt.subplots(1,1, sharex =False, sharey=False, figsize=(10,8))
plot2 = sns.scatterplot("I", "Ipref", data=out, ax=axarr)
fig.savefig("nprefTestI.png")
