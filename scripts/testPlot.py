import pandas
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import sys
import mating
from createPop import *
import os

newFile = str(len(os.listdir("output/freqs"))+1)

with open("output/paramValues.tsv", "r") as paramFile:
    nameList = paramFile.readline()
    nameList = nameList.strip()
    nameList = nameList.split("\t")

with open("output/paramValues.tsv", "a") as paramFile:
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
    data = mating.runSim("long")
    with open("output/raw/sim{}.tsv".format(newFile), "w") as dataFile:
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

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["A", "I", "O"], palette=["blue", "green", "red"], data=table, legend=False)
sns.despine()
plot.set(ylim=(-0.01,1.01))
plot.fig.suptitle("Phenotype Frequencies", y=1)
plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
plot.set(xlabel="Generation", ylabel="Frequency")
plot.savefig("output/freqs/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["M", "F", "T"], palette =["blue", "red", "black"], data=table, legend=False)
sns.despine()
plot.fig.suptitle("Population Size", y=1)
plt.legend(labels=("Male", "Female", "Total"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Sex")
plot.set(xlabel="Generation", ylabel="Population")
plot.savefig("output/popSize/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["MalF", "FemF"], palette=["blue", "red"], data=table, legend=False)
sns.despine()
plot.set(ylim=(-0.01,1.01))
plot.fig.suptitle("Fertility", y=1)
plt.legend(labels=("Male", "Female"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Sex")
plot.set(xlabel="Generation", ylabel="fertility")
plot.savefig("output/fertility/sim{}S.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order=["AF", "IF", "OF"], palette=["blue", "green", "red"], data=table,legend=False)
sns.despine()
plot.set(ylim=(-0.01,1.01))
plot.fig.suptitle("Fertility", y=1)
plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
plot.set(xlabel="Generation", ylabel="fertility")
plot.savefig("output/fertility/sim{}P.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["APref", "IPref", "OPref"], palette = ["blue", "green", "red"], data=table, legend=False)
sns.despine()
plot.set(ylim=(0,1.01))
plot.fig.suptitle("Phenotype Preference Post-Mating", y=1)
plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
plot.set(xlabel="Generation", ylabel="Preference")
plot.savefig("output/prefs/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["Contacts", "Matings", "MMContacts"], palette = ["red", "green", "blue"], data=table, legend=False)
sns.despine()
plot.fig.suptitle("Interactions", y=1)
plt.legend(labels=("M-F Interactions", "Copulations", "M-M Interactions"), loc="upper left", bbox_to_anchor=(1.05, 0.5))
plot.set(xlabel="Generation", ylabel="Number of Interactions")
plot.savefig("output/contact/sim{}.png".format(newFile), bbox_inches="tight")

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["preAPref", "preIPref", "preOPref"], palette = ["blue", "green", "red"], data=table, legend=False)
sns.despine()
plot.set(ylim=(0,1.01))
plot.fig.suptitle("Phenotype Frequency Pre-mating", y=1)
plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
plot.set(xlabel="Generation", ylabel="Preference")
plot.savefig("output/prePrefs/sim{}.png".format(newFile), bbox_inches="tight")
