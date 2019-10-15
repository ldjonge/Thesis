import numpy as np
import pandas as pd
import sys
from standardFunctions import wavg

def reshape(data):
    parameters = data.Pheno.unique()
    Runs = data.Run.unique()
    Gens = data.Gen.unique()
    added = []

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
    data[["M", "F", "T", "Matings", "Contacts"]] = data[["M", "F", "T", "Matings", "Contacts"]].astype('Int64')
    totalPopDict = {}
    for col in data.select_dtypes('Int64').columns:
        totalPopDict[col] = data.groupby(["Run", "Gen"])[col].sum()

    for col in ["A", "I", "O", "FFec"]:
        totalPopDict[col] = data.groupby(["Run", "Gen"]).apply(wavg, col, "F")

    for col in ["preAPref", "preIPref", "preOPref", "APref", "IPref", "OPref"]:
        totalPopDict[col] = data.groupby(["Run", "Gen"]).apply(wavg, col, "M")

    #totalPopDict["Run", "Gen"] = data.groupby(["Run", "Gen"]).groups.keys()
    totalPopDict["Pop"] = ["Total" for i in range(len(totalPopDict["A"]))]
    totalPop = pd.DataFrame(totalPopDict)
    totalPop = totalPop.reset_index()
    data = data.reset_index()

    data = data.append(totalPop, sort=False)

    data['matingSuccess'] = (data['Matings']/data['Contacts']).astype(float)
    #data['misIdent'] = (data['MMContacts']/(data['Contacts']+data['MMContacts'])).astype(float)

    data.loc[data['F']==0, 'A'] = np.nan
    data.loc[data['F']==0, 'I'] = np.nan
    data.loc[data['F']==0, 'O'] = np.nan

    return data


if __name__=="__main__":

    if len(sys.argv)>1:
        data = pd.read_csv(sys.argv[1], sep="\t")
    else:
        data = pd.read_csv("multiOutput/summary/summary.tsv", sep="\t")

    data = reshape(data)

    data.to_csv("multiOutput/summary/newData1.tsv", sep="\t", index=False)
