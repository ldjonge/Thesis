import pandas
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import sys
import mating
import newMultiPop
from createPop import *
import os

newFile = str(len(os.listdir("multiOutput/freqs"))+1)
print("Simulation ", newFile)

with open("multiOutput/paramValues.tsv", "r") as paramFile:
    nameList = paramFile.readline()
    nameList = nameList.strip()
    nameList = nameList.split("\t")

with open("multiOutput/paramValues.tsv", "a") as paramFile:
    paramValues = readParams()
    valueList = [newFile]
    for i in nameList[1:]:
        try:
            valueList.append(str(paramValues[i]))
        except KeyError:
            print(i)
            valueList.append(" ")
    valueList = "\t".join(valueList)
    print(valueList, file=paramFile)

if len(sys.argv) > 1:
    table = pandas.read_csv(sys.argv[1], sep="\t")

else:
    data = newMultiPop.runSim("long")
    with open("multiOutput/raw/sim{}.tsv".format(newFile), "w") as dataFile:
        for line in data:
            new = []
            for i in line:
                new.append(str(i))
            print("\t".join(new), file=dataFile)
    headers = data.pop(0)
    table = pandas.DataFrame(data, columns=headers)
#print(table)


plt.figure(figsize=(10,8))
sns.set_style("white")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["A", "I", "O"], style="Pop", palette=["blue", "green", "red"], data=table, legend='brief')
sns.despine()
plot.set(ylim=(-0.01,1.01))
plot.fig.suptitle("Phenotype Frequencies", y=1)
#plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
plot.set(xlabel="Generation", ylabel="Frequency")
plot.savefig("multiOutput/freqs/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["M", "F", "T"], style="Pop", palette =["blue", "red", "black"], data=table, legend='brief')
sns.despine()
plot.fig.suptitle("Population Size", y=1)
#plt.legend(labels=("Male", "Female", "Total"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Sex")
plot.set(xlabel="Generation", ylabel="Population")
plot.savefig("multiOutput/popSize/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["MalF", "FemF"], style="Pop", palette=["blue", "red"], data=table, legend='brief')
sns.despine()
plot.set(ylim=(-0.01,1.01))
plot.fig.suptitle("Fecundity", y=1)
#plt.legend(labels=("Male", "Female"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Sex")
plot.set(xlabel="Generation", ylabel="Fecundity")
plot.savefig("multiOutput/fecundity/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["APref", "IPref", "OPref"], style="Pop", palette = ["blue", "green", "red"], data=table, legend='brief')
sns.despine()
plot.set(ylim=(0,1.01))
plot.fig.suptitle("Phenotype Preference Post-Mating", y=1)
# plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
plot.set(xlabel="Generation", ylabel="Preference")
plot.savefig("multiOutput/prefs/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["Contacts", "Matings", "MMContacts"], style="Pop", palette = ["red", "green", "blue"], data=table, legend='brief')
sns.despine()
plot.fig.suptitle("Interactions", y=1)
# plt.legend(labels=("M-F Interactions", "Copulations", "M-M Interactions"), loc="upper left", bbox_to_anchor=(1.05, 0.5))
plot.set(xlabel="Generation", ylabel="Number of Interactions")
plot.savefig("multiOutput/contact/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["preAPref", "preIPref", "preOPref"], style="Pop", palette = ["blue", "green", "red"], data=table, legend='brief')
sns.despine()
plot.set(ylim=(0,1.01))
plot.fig.suptitle("Phenotype Frequency Pre-mating", y=1)
# plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
plot.set(xlabel="Generation", ylabel="Preference")
plot.savefig("multiOutput/prePrefs/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["migrations"], style="Pop", palette = ["green"], data=table, legend="brief")
sns.despine()
plot.fig.suptitle("Migration", y=1)
plot.set(xlabel="Generation")
plot.savefig("multiOutput/migrations/sim{}.png".format(newFile), bbox_inches="tight")
