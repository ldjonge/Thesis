import pandas
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) > 1 and sys.argv[1].lower() == "multi":
    multi = True
    data = pandas.read_csv("multiOutput/summary/summary.tsv", sep="\t")
else:
    multi = False
    data = pandas.read_csv("newOutputAgain/summary/summary.tsv", sep="\t")
parameters = data.Pheno.unique()


dataList = []
if multi:
    for par in parameters:
        parDF = data.loc[data['Pheno']==par]
        parDF = parDF[['Pop', 'Run', 'Gen', 'Value']]
        parDF = parDF.groupby(['Pop', 'Gen', 'Run']).Value.sum().to_frame()
        parDF.rename(index=str, columns={"Value": par}, inplace=True)
        dataList.append(parDF)

    dataReshape = dataList[0].join(dataList[1], on=["Pop", "Gen", "Run"])
    for df in dataList[2:]:
        dataReshape = dataReshape.join(df, on=["Pop", "Gen", "Run"])
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
print(data['matingSuccess'].min(), data['matingSuccess'].max(), data['matingSuccess'].mean())
data['misIdent'] = (data['MMContacts']/(data['Contacts']+data['MMContacts'])).astype(float)

print(data[data['matingSuccess']==data['matingSuccess'].max()])

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
    fig.savefig("multiOutput/summary/correlation.png")
else:
    fig.savefig("newOutputAgain/summary/correlation.png")
#print(data.head())
sns.set_style("white")

#Fertility vs Androchrome Frequency
fig, axarr = plt.subplots(1,2, sharey=True, figsize=(10,8))

if multi:
    plot1 = sns.regplot(x="A", y="MalF", style="Pop", data=data,ax=axarr[0])
    plot1.set_title("Male")
    plot2 = sns.regplot(x="A", y="FemF", style="Pop", data=data, ax=axarr[1])
    plot2.set_title("Female")
    sns.despine()
    fig.suptitle("Fecundity")
    fig.savefig("multiOutput/summary/sexFec.png")
else:
    plot1 = sns.regplot(x="A", y="MalF", data=data,ax=axarr[0])
    plot1.set_title("Male")
    plot2 = sns.regplot(x="A", y="FemF", data=data, ax=axarr[1])
    plot2.set_title("Female")
    sns.despine()
    fig.suptitle("Fecundity")
    fig.savefig("newOutputAgain/summary/sexFec.png")

#Fertility vs Frequency
if not multi:
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
if multi:
    plot1 = sns.regplot(x="A", y="preAPref", style="Pop", data=data, ax=axarr[0,0])
    plot1.set(ylabel="")
    plot1.set_title("Pre-Mating")
    plot2 = sns.regplot(x="A", y="APref", style="Pop", data=data, ax=axarr[0,1])
    plot2.set(ylabel="")
    plot2.set_title("Post-Mating")
    plot3 = sns.regplot(x="I", y="preIPref", style="Pop", data=data, ax=axarr[1,0])
    plot3.set(ylabel="")
    plot4 = sns.regplot(x="I", y="IPref", style="Pop", data=data, ax=axarr[1,1])
    plot4.set(ylabel="")
    plot5 = sns.regplot(x="O", y="preOPref", style="Pop", data=data, ax=axarr[2,0])
    plot5.set(ylabel="")
    plot6 = sns.regplot(x="O", y="OPref", style="Pop", data=data, ax=axarr[2,1])
    plot6.set(ylabel="")
    sns.despine()

    fig.savefig("multiOutput/summary/prefFreq.png")
else:
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
if multi:
    plot = sns.regplot(x="A", y="misIdent", style="Pop", data=data)
else:
    plot = sns.regplot(x="A", y="misIdent", data=data)
sns.despine()
if multi:
    plt.savefig("multiOutput/summary/maleMating.png")
else:
    plt.savefig("newOutputAgain/summary/maleMating.png")
#plot = sns.lmplot(x="A", y="preAPref", data=data, legend=False)

plt.figure(figsize=(10,8))
if multi:
    plot = sns.regplot(x="M", y="F", style="Pop", data=data)
else:
    plot = sns.regplot(x="M", y="F", data=data)

sns.despine()
if multi:
    plt.savefig("multiOutput/summary/popDist.png")
else:
    plt.savefig("newOutputAgain/summary/popDist.png")
