import pandas
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import sys
#import mating
#import newMultiPop
import newLearning as nL
#from createPop import *
import os
from readFiles import *

if len(sys.argv) > 1:
    table = pandas.read_csv(sys.argv[1], sep="\t")
    newFile= str(123)

else:
    newFile = str(len(os.listdir("multiOutput/raw"))+1)
    print("Simulation ", newFile)
    data = nL.runSim()
    with open("multiOutput/raw/sim{}.tsv".format(newFile), "w") as dataFile:
        for line in data:
            new = []
            for i in line:
                new.append(str(i))
            print("\t".join(new), file=dataFile)
    headers = data.pop(0)
    table = pandas.DataFrame(data, columns=headers)
#print(table)

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
            #print(i)
            valueList.append(" ")
    valueList = "\t".join(valueList)
    print(valueList, file=paramFile)

plt.figure(figsize=(10,8))
sns.set_style("white")
for pop in table.Pop.unique():
    print(pop)
    plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["A", "I", "O"], palette=["blue", "green", "red"], data=table[table["Pop"]==pop], legend='brief')
    sns.despine()
    plot.set(ylim=(-0.01,1.01))
    plot.fig.suptitle("Phenotype Frequencies", y=1)
    #plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
    plot.set(xlabel="Generation", ylabel="Frequency")
    plot.savefig("multiOutput/freqs/sim{}pop{}.png".format(newFile, pop), bbox_inches="tight")
    plt.close()

    plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["T"], palette =["black"], data=table[table["Pop"]==pop], legend='brief')
    sns.despine()
    plot.fig.suptitle("Population Size", y=1)
    plot.set(ylim=(0, max(table[table["Pheno"]=="T"]["Value"])))
    #plt.legend(labels=("Male", "Female", "Total"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Sex")
    plot.set(xlabel="Generation", ylabel="Population")
    plot.savefig("multiOutput/popSize/sim{}pop{}.png".format(newFile, pop), bbox_inches="tight")
    plt.close()

    plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["AFec", "IFec", "OFec"], palette=["blue", "green", "red"], data=table[table["Pop"]==pop], legend='brief')
    sns.despine()
    plot.set(ylim=(0,max(table[table["Pheno"]=="FFec"]["Value"])))
    plot.fig.suptitle("Fecundity", y=1)
    #plt.legend(labels=("Male", "Female"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Sex")
    plot.set(xlabel="Generation", ylabel="Fecundity")
    plot.savefig("multiOutput/fecundity/sim{}pop{}.png".format(newFile, pop), bbox_inches="tight")
    plt.close()

    plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["APref", "IPref", "OPref"], palette = ["blue", "green", "red"], data=table[table["Pop"]==pop], legend='brief')
    sns.despine()
    plot.set(ylim=(-0.01,1.01))
    plot.fig.suptitle("Phenotype Preference Post-Mating", y=1)
    # plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
    plot.set(xlabel="Generation", ylabel="Preference")
    plot.savefig("multiOutput/prefs/sim{}pop{}.png".format(newFile, pop), bbox_inches="tight")
    plt.close()

    plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["Contacts", "Matings"], palette = ["red", "green"], data=table[table["Pop"]==pop], legend='brief')
    sns.despine()
    plot.fig.suptitle("Interactions", y=1)
    plot.set(ylim=(0, max(table[table["Pheno"]=="Contacts"]["Value"])))
    # plt.legend(labels=("M-F Interactions", "Copulations", "M-M Interactions"), loc="upper left", bbox_to_anchor=(1.05, 0.5))
    plot.set(xlabel="Generation", ylabel="Number of Interactions")
    plot.savefig("multiOutput/contact/sim{}pop{}.png".format(newFile, pop), bbox_inches="tight")
    plt.close()

    plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["preAPref", "preIPref", "preOPref"], palette = ["blue", "green", "red"], data=table[table["Pop"]==pop], legend='brief')
    sns.despine()
    plot.set(ylim=(0,1.01))
    plot.fig.suptitle("Phenotype Frequency Pre-mating", y=1)
    # plt.legend(labels=("Androchrome", "Infuscans", "Infuscans-Obsoleta"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Phenotype")
    plot.set(xlabel="Generation", ylabel="Preference")
    plot.savefig("multiOutput/prePrefs/sim{}pop{}.png".format(newFile, pop), bbox_inches="tight")
    plt.close()

    plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["MMSucc"], palette =["black"], data=table[table["Pop"]==pop], legend='brief')
    sns.despine()
    plot.fig.suptitle("Mating Success", y=1)
    #plt.legend(labels=("Male", "Female", "Total"), loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Sex")
    plot.set(xlabel="Generation", ylabel="Population")
    plot.savefig("multiOutput/mSucc/sim{}pop{}.png".format(newFile, pop), bbox_inches="tight")
    plt.close()

    """
    plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["Migrations"], palette = ["green"], data=table[table["Pop"]==pop], legend="brief")
    sns.despine()
    plot.fig.suptitle("Migration", y=1)
    plot.set(xlabel="Generation")
    plot.savefig("multiOutput/migrations/sim{}pop{}.png".format(newFile, pop), bbox_inches="tight")
    plt.close()

    plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["Fertilised Migrations"], palette = ["green"], data=table[table["Pop"]==pop], legend="brief")
    sns.despine()
    plot.fig.suptitle("Migration", y=1)
    plot.set(xlabel="Generation")
    plot.savefig("multiOutput/fertmigrations/sim{}pop{}.png".format(newFile, pop), bbox_inches="tight")
    plt.close()
    """
