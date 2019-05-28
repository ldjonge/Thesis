import pandas
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import numpy as np
import sys
from heatmap import *

if len(sys.argv) > 1 and sys.argv[1].lower() == "multi":
    multi = True
    data = pandas.read_csv("multiOutput/summary/summary/summary.tsv", sep="\t")
else:
    multi = False
    data = pandas.read_csv("newOutputAgain/summary/summary.tsv", sep="\t")
parameters = data.Pheno.unique()


dataList = []
if multi:
    for par in parameters:
        parDF = data.loc[data['Pheno']==par]
        parDF = parDF[['Pop', 'Run', 'Gen', 'Value']]
        parDF = parDF.groupby(['Pop', 'Run', 'Gen']).Value.sum().to_frame()
        parDF.rename(index=str, columns={"Value": par}, inplace=True)
        dataList.append(parDF)

    dataReshape = dataList[0].join(dataList[1], on=["Pop", "Run", "Gen"])
    for df in dataList[2:]:
        dataReshape = dataReshape.join(df, on=["Pop", "Run", "Gen"])

else:
    for par in parameters:
        parDF = data.loc[data['Pheno']==par]
        parDF = parDF[['Run', 'Gen', 'Value']]
        parDF = parDF.groupby(['Gen', 'Run']).Value.sum().to_frame()
        parDF.rename(index=str, columns={"Value": par}, inplace=True)
        dataList.append(parDF)

    dataReshape = dataList[0].join(dataList[1], on=["Gen", "Run"])
    for df in dataList[2:]:
        dataReshape = dataReshape.join(df, on=["Gen", "Run"])

#print(dataReshape.head())
data = dataReshape

data['matingSuccess'] = (data['Matings']/data['Contacts']).astype(float)
#print(data['matingSuccess'].min(), data['matingSuccess'].max(), data['matingSuccess'].mean())
data['misIdent'] = (data['MMContacts']/(data['Contacts']+data['MMContacts'])).astype(float)

#print(data[data['matingSuccess']==data['matingSuccess'].max()])

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
if multi:
    fig.savefig("multiOutput/summary/summary/correlation.png")
    plt.figure(figsize=(10,10))
    corrplot(data.corr())
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.savefig("multiOutput/summary/summary/corrHeatMap.png")
else:
    fig.savefig("newOutputAgain/summary/correlation.png")
    plt.figure(figsize=(10,10))
    corrplot(data.corr())
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.savefig("newOutputAgain/summary/corrHeatMap.png")
#print(data.head())


#Fertility vs Androchrome Frequency

if multi:
    pops = dataReshape.groupby("Pop")
    for pop, popData in pops:
        print(pop)
        print(popData.describe())
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
        plt.close(fig='all')
else:
    sns.set_style("white")
    fig, axarr = plt.subplots(1,2, sharey=True, figsize=(10,8))
    plot1 = sns.regplot(x="A", y="MalF", data=data,ax=axarr[0])
    plot1.set_title("Male")
    plot2 = sns.regplot(x="A", y="FemF", data=data, ax=axarr[1])
    plot2.set_title("Female")
    sns.despine()
    fig.suptitle("Fecundity")
    fig.savefig("newOutputAgain/summary/sexFec.png")

    #Fertility vs Frequency
    fig, axarr = plt.subplots(1,3, sharey=True, figsize=(10,8))
    plot1 = sns.regplot(x="A", y="AF", data=data, ax=axarr[0])
    plot1.set_title("A")
    plot2=sns.regplot(x="I", y="IF", data=data, ax=axarr[1])
    plot2.set_title("I")
    plot3=sns.regplot(x="O", y="OF", data=data, ax=axarr[2])
    plot3.set_title("O")
    sns.despine()
    fig.suptitle("Fecundity")
    fig.savefig("newOutputAgain/summary/phenFec.png")

#Preference vs Frequency
    fig, axarr = plt.subplots(3,2, sharex ='col', sharey=False, figsize=(10,8))
    fig.suptitle("Preference")
    plot1 = sns.regplot(x="A", y="preAPref", data=data, ax=axarr[0,0])
    plot1.set(ylabel="")
    plot1.set_title("Pre-Mating")
    plot2 = sns.regplot(x="A", y="APref", data=data, ax=axarr[0,1])
    plot2.set(ylabel="")
    plot2.set_title("Post-Mating")
    plot3 = sns.regplot(x="I", y="preIPref", data=data, ax=axarr[1,0])
    plot3.set(ylabel="")
    plot4 = sns.regplot(x="I", y="IPref", data=data, ax=axarr[1,1])
    plot4.set(ylabel="")
    plot5 = sns.regplot(x="O", y="preOPref", data=data, ax=axarr[2,0])
    plot5.set(ylabel="")
    plot6 = sns.regplot(x="O", y="OPref", data=data, ax=axarr[2,1])
    plot6.set(ylabel="")
    sns.despine()
    fig.savefig("newOutputAgain/summary/prefFreq.png")


    plt.figure(figsize=(10,8))
    plot = sns.regplot(x="A", y="misIdent", data=data)
    sns.despine()
    plt.savefig("newOutputAgain/summary/maleMating.png")


    plt.figure(figsize=(10,8))
    plot = sns.regplot(x="M", y="F", data=data)

    sns.despine()

    plt.savefig("newOutputAgain/summary/popDist.png")
