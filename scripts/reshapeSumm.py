import numpy as np
import pandas as pd
import sys


if len(sys.argv)>1:
    data = pd.read_csv(sys.argv[1], sep="\t")
else:
    data = pd.read_csv("multiOutput/summary/summary/summary.tsv", sep="\t")

parameters = data.Pheno.unique()
Runs = data.Run.unique()
Gens = data.Gen.unique()
added = []

def wavg(group, avgName, weightName):
    d = group[avgName]
    w = group[weightName]
    try:
        return (d*w).sum() / w.sum()
    except ZeroDivisionError:
        return d.mean()

""" Figure out how to cleanly apply weighted averages for the populations where needed!! """
"""
for par in parameters:
    parData = data.loc[data['Pheno']==par]
    print(par)
    for run in Runs:
        runData = parData.loc[parData['Run']==run]
        for gen in Gens:
            genData = runData.loc[runData['Gen']==gen]
            if par in ["M", "F", "T", "migrations", "Matings", "Contacts", "MMContacts"]:
                added.append(pd.DataFrame({"Pop":["Total"], "Gen":[gen], "Run":[run], "Pheno": [par], "Value": [genData.Value.sum()]}))
            else:
                F =
                added.append(pd.DataFrame({"Pop":["Total"], "Gen":[gen], "Run":[run], "Pheno": [par], "Value": [genData.Value.mean()]}))
data = data.append(added, ignore_index=True, sort=True)
"""
dataList=[]
for par in parameters:
    parDF = data.loc[data['Pheno']==par]
    parDF = parDF[['Pop', 'Run', 'Gen', 'Value']]
    parDF = parDF.groupby(['Pop', 'Run', 'Gen']).Value.sum().to_frame()
    parDF.rename(index=str, columns={"Value": par}, inplace=True)
    dataList.append(parDF)

dataReshape = dataList[0].join(dataList[1], on=["Pop", "Run", "Gen"])
for df in dataList[2:]:
    dataReshape = dataReshape.join(df, on=["Pop", "Run", "Gen"])
data = dataReshape

data[["M", "F", "T", "migrations", "Matings", "Contacts", "MMContacts"]] = data[["M", "F", "T", "migrations", "Matings", "Contacts", "MMContacts"]].astype('Int64')
totalPopDict = {}
for col in data.select_dtypes('Int64').columns:
    totalPopDict[col] = data.groupby(["Run", "Gen"])[col].sum()

for col in ["A", "I", "O", "FemF"]:
    totalPopDict[col] = data.groupby(["Run", "Gen"]).apply(wavg, col, "F")

for col in ["preAPref", "preIPref", "preOPref", "APref", "IPref", "OPref", "MalF"]:
    totalPopDict[col] = data.groupby(["Run", "Gen"]).apply(wavg, col, "M")

#totalPopDict["Run", "Gen"] = data.groupby(["Run", "Gen"]).groups.keys()
totalPopDict["Pop"] = ["Total" for i in range(len(totalPopDict["A"]))]
totalPop = pd.DataFrame(totalPopDict)
totalPop = totalPop.reset_index()
data = data.reset_index()
data = data.append(totalPop, sort=False)

data['matingSuccess'] = (data['Matings']/data['Contacts']).astype(float)
data['misIdent'] = (data['MMContacts']/(data['Contacts']+data['MMContacts'])).astype(float)

data['M'].replace(0, np.nan, inplace=True)
data['F'].replace(0, np.nan, inplace=True)
data['T'].replace(0, np.nan, inplace=True)
#print(data.describe())


data.to_csv("multiOutput/summary/summary/newData.tsv", sep="\t", index=False)
