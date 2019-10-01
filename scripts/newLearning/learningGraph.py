from pops import *
from baseFunctions import *
from readFiles import *
from record import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

class simpMale:
    def __init__(self):
        self.prefs = {"A": 1, "I": 1, "O": 1}

    def learning(self, pop, prefEff):
        chance = [fem.vis for fem in pop]
        chSum = sum(chance)
        chance = [i/chSum for i in chance]
        for fem in np.random.choice(pop, size=min(len(pop), 10), p=chance, replace=False):
            # if random.random() < fem.mSucc:
            self.prefs[fem.phenotype] *= prefEff
            #else:
            #    self.prefs[fem.phenotype] *= params["failPrefEff"]
        valsum = sum(self.prefs.values())
        for val in self.prefs.values():
            val = val/valsum

freqOptions = []
for i in range(101):
    for j in range(101-i):
        k = 100-i-j
        freqOptions.append((i,j,k))

outTable = [["A", "I", "O", "Apref", "Ipref", "Opref"]]
prefEff = 2
mimicEff = "variable2"
for freqs in freqOptions:
    malePop = [simpMale() for i in range(100)]
    femPop = [Female("p", "p") for i in range(freqs[0])]
    for i in range(freqs[1]):
        femPop.append(Female("q", "q"))
    for i in range(freqs[2]):
        femPop.append(Female("r", "r"))
    for female in femPop:
        if female.phenotype == "A":
            female.vis = (freqs[0]/(100+freqs[0]))**0.25
        else:
            female.vis = 1
    for male in malePop:
        male.learning(femPop, prefEff)
    prefs = recordPref([malePop, femPop])
    outTable.append(list(freqs)+prefs)

outTable = pd.DataFrame(outTable[1:], columns=outTable[0])
#print(outTable.head())

plt.figure(figsize=(10,8))
sns.set_style("white")

fig, axarr = plt.subplots(1,1, figsize=(20,16))
plot = sns.scatterplot(x="A", y="Apref", data=outTable, ax=axarr)
fig.suptitle("A Preference, Strength = {}, Mimicry={}".format(prefEff, mimicEff))
plot.set(ylim=(0,1))
plot.set(xlabel="Frequency", ylabel="Preference")
plt.savefig("learningCurves/APref_{}_mimic_{}_highDens.png".format(prefEff, mimicEff), bbox_inches="tight")

fig, axarr = plt.subplots(1,1, figsize=(20,16))
plot = sns.scatterplot(x="I", y="Ipref", data=outTable, ax=axarr)
fig.suptitle("I Preference, Strength = {}, Mimicry={}".format(prefEff, mimicEff))
plot.set(ylim=(0,1))
plot.set(xlabel="Frequency", ylabel="Preference")
plt.savefig("learningCurves/IPref_{}_mimic_{}_highDens.png".format(prefEff, mimicEff), bbox_inches="tight")

fig, axarr = plt.subplots(1,1, figsize=(20,16))
plot = sns.scatterplot(x="O", y="Opref", data=outTable, ax=axarr)
fig.suptitle("O Preference, Strength = {}, Mimicry={}".format(prefEff, mimicEff))
plot.set(ylim=(0,1))
plot.set(xlabel="Frequency", ylabel="Preference")
plt.savefig("learningCurves/OPref_{}_mimic_{}_highDens.png".format(prefEff, mimicEff), bbox_inches="tight")
plt.close()
