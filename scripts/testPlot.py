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

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["A", "I", "O"], palette=["blue", "green", "red"], data=table)
sns.despine()
plot.set(ylim=(None,1.01))
plot.savefig("output/freqs/sim{}.png".format(newFile))

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["M", "F", "T"], palette =["blue", "red", "black"], data=table)
sns.despine()
plot.savefig("output/popSize/sim{}.png".format(newFile))

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["MalF", "FemF"], palette=["blue", "red"], data=table)
sns.despine()
plot.set(ylim=(None,1.01))
plot.savefig("output/fecundity/sim{}.png".format(newFile))

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["APref", "IPref", "OPref"], palette = ["blue", "green", "red"], data=table)
sns.despine()
plot.set(ylim=(None,1.01))
plot.savefig("output/prefs/sim{}.png".format(newFile))

plot = sns.relplot(x="Gen", y="Value", kind="line", hue="Pheno", hue_order = ["Contacts", "Matings", "MMContacts"], palette = ["red", "green", "blue"], data=table)
sns.despine()
plot.savefig("output/contact/sim{}.png".format(newFile))
