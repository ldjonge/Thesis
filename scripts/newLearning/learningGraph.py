from pops import *
from baseFunctions import *
from readFiles import *
from record import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

class simpMale:
    def __init__(self):
        self.prefs = {"A": 0.02, "I": 0.02, "O": 0.02}

    def learning(self, pop, prefEff):
        for fem in np.random.choice(pop, size=min(len(pop), 20), replace=True):
            if type(fem) == simpMale:
                self.prefs["A"] *= 0.9
            elif type(fem) == Female:
                if random.random() < 0.4:
                    self.prefs[fem.phenotype] *= prefEff
                else:
                    self.prefs[fem.phenotype] *= 0.9 #params["failPrefEff"]
        valsum = sum(self.prefs.values())+3
        for key in self.prefs.keys():
            self.prefs[key] = 0.7*(((1+self.prefs[key]))/valsum)+0.1

freqOptions = []
for i in range(101):
    for j in range(101-i):
        k = 100-i-j
        freqOptions.append((i,j,k))

outTable = [["A", "I", "O", "Apref", "Ipref", "Opref"]]
prefEff = 5
mimicEff = 1 #"Mvariable"
for freqs in freqOptions:
    prefs = [0, 0, 0]
    malePop = [simpMale() for i in range(50)]
    femPop = [Female("p", "p") for i in range(freqs[0])]
    for i in range(freqs[1]):
        femPop.append(Female("q", "q"))
    for i in range(freqs[2]):
        femPop.append(Female("r", "r"))
    for male in malePop:
        male.learning(femPop+malePop, prefEff)
        prefs[0] += male.prefs["A"]
        prefs[1] += male.prefs["I"]
        prefs[2] += male.prefs["O"]
    prefs = [i/len(malePop) for i in prefs]
    outTable.append(list(freqs)+prefs)

outTable = pd.DataFrame(outTable[1:], columns=outTable[0])
outTable["AProb"]=outTable["A"]*outTable["Apref"]/100
print(outTable.head())

plt.figure(figsize=(10,8))
sns.set_style("white")
sns.set_context("talk")

fig, axarr = plt.subplots(1,1, figsize=(10,8))
plot = sns.scatterplot(x="A", y="Apref", data=outTable, ax=axarr)
#fig.suptitle("A Preference, Strength = {}, Mimicry={}".format(prefEff, mimicEff))
plot.set(ylim=(0,1))
plot.set(xlabel="Frequency", ylabel="Preference")
print("hello World")
plt.savefig("newLearningCurves/newAPref_{}.png".format(prefEff), bbox_inches="tight")

fig, axarr = plt.subplots(1,1, figsize=(10,8))
plot = sns.scatterplot(x="I", y="Ipref", data=outTable, ax=axarr)
#fig.suptitle("I Preference, Strength = {}, Mimicry={}".format(prefEff, mimicEff))
plot.set(ylim=(0,1))
plot.set(xlabel="Frequency", ylabel="Preference")
plt.savefig("newLearningCurves/newIPref_{}".format(prefEff), bbox_inches="tight")

fig, axarr = plt.subplots(1,1, figsize=(10,8))
plot = sns.scatterplot(x="O", y="Opref", data=outTable, ax=axarr)
#fig.suptitle("O Preference, Strength = {}, Mimicry={}".format(prefEff, mimicEff))
plot.set(ylim=(0,1))
plot.set(xlabel="Frequency", ylabel="Preference")
plt.savefig("newLearningCurves/newOPref_{}.png".format(prefEff), bbox_inches="tight")
plt.close()

fig, axarr = plt.subplots(1,1, figsize=(20,16))
plot = sns.scatterplot(x="A", y="AProb", data=outTable, ax=axarr)
plot.set(ylim=(0,1))
plt.savefig("newLearningCurves/AProb{}.png".format(prefEff), bbox_inches="tight")
