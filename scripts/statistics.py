import numpy as np
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from heatmap import *
import statsmodels.api as sm
import statsmodels.formula.api as smf
import sys

if len(sys.argv)>1:
    data = pd.read_csv(sys.argv[1], sep="\t", dtype={'Pop':'str', 'Run':'int', 'Gen':'int', 'preAPref':'float', 'preIPref':'float', 'preOPref':'float', 'Matings':'float','Contacts':'float', 'MMContacts':'float', 'MalF':'float', 'FemF':'float', 'APref':'float', 'IPref':'float', 'OPref':'float','migrations':'float', 'A':'float', 'I':'float', 'O':'float', 'M':'float', 'F':'float', 'T':'float', 'matingSuccess':'float', 'misIdent':'float'})
else:
    data = pd.read_csv("multiOutput/summary/summary/newData.tsv", sep="\t",  dtype={'Pop':'str', 'Run':'int', 'Gen':'int', 'preAPref':'float', 'preIPref':'float', 'preOPref':'float', 'Matings':'float','Contacts':'float', 'MMContacts':'float', 'MalF':'float', 'FemF':'float', 'APref':'float', 'IPref':'float', 'OPref':'float','migrations':'float', 'A':'float', 'I':'float', 'O':'float', 'M':'float', 'F':'float', 'T':'float', 'matingSuccess':'float', 'misIdent':'float'})

data.set_index(["Pop", "Run", "Gen"], inplace=True)

pops = data.groupby("Pop")



""" PLOTS """
sns.set_context("talk")

corr = data.corr()
fig = plt.figure(figsize=(15,12))
fig.suptitle("Correlation Matrix")
ax = fig.add_subplot(111)
cax = ax.matshow(corr,cmap='coolwarm', vmin=-1, vmax=1)
fig.colorbar(cax)
ticks = np.arange(0,len(data.columns),1)
ax.set_xticks(ticks)
plt.xticks(rotation=90)
ax.set_yticks(ticks)
ax.set_xticklabels(data.columns)
ax.set_yticklabels(data.columns)

fig.savefig("multiOutput/summary/summary/correlation.png")
plt.figure(figsize=(10,10))
corrplot(data.corr())
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.savefig("multiOutput/summary/summary/corrHeatMap.png")

for pop, popData in pops:
    sns.set_style("darkgrid")
    plt.figure(figsize=(10,10))
    corrplot(popData.corr())
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.savefig("multiOutput/summary/summary/corrHeatMap{}.png".format(pop))
    sns.set_style("white")
    fig, axarr = plt.subplots(1,2, sharey=True, figsize=(10,8))
    plot1 = sns.regplot(x="A", y="MalF", data=popData,ax=axarr[0])
    plot1.set_title("Male")
    plot2 = sns.regplot(x="A", y="FemF", data=popData, ax=axarr[1])
    plot2.set_title("Female")
    sns.despine()
    fig.suptitle("Fecundity")
    fig.savefig("multiOutput/summary/summary/sexFec{}.png".format(pop))

    #Preference vs Frequency
    fig, axarr = plt.subplots(3,2, sharex ='col', sharey=False, figsize=(10,8))
    fig.suptitle("Preference")
    plot1 = sns.regplot(x="A", y="preAPref", data=popData, ax=axarr[0,0])
    plot1.set(ylabel="")
    plot1.set_title("Pre-Mating")
    plot2 = sns.regplot(x="A", y="APref", data=popData, ax=axarr[0,1])
    plot2.set(ylabel="")
    plot2.set_title("Post-Mating")
    plot3 = sns.regplot(x="I", y="preIPref", data=popData, ax=axarr[1,0])
    plot3.set(ylabel="")
    plot4 = sns.regplot(x="I", y="IPref", data=popData, ax=axarr[1,1])
    plot4.set(ylabel="")
    plot5 = sns.regplot(x="O", y="preOPref", data=popData, ax=axarr[2,0])
    plot5.set(ylabel="")
    plot6 = sns.regplot(x="O", y="OPref", data=popData, ax=axarr[2,1])
    plot6.set(ylabel="")
    sns.despine()
    fig.savefig("multiOutput/summary/summary/prefFreq{}.png".format(pop))

    plt.figure(figsize=(10,8))
    plot = sns.regplot(x="A", y="misIdent", data=popData)
    sns.despine()
    plt.savefig("multiOutput/summary/summary/maleMating{}.png".format(pop))

    plt.figure(figsize=(10,8))
    plot = sns.regplot(x="M", y="F", data=popData)
    sns.despine()
    plt.savefig("multiOutput/summary/summary/popDist{}.png".format(pop))


    plt.figure(figsize=(10,8))
    plot = sns.regplot(x="T", y="migrations", data=popData)
    sns.despine()
    plt.savefig("multiOutput/summary/summary/migr{}.png".format(pop))
    plt.close(fig='all')
