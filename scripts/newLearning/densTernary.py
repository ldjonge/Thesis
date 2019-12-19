import numpy as np
import pandas as pd
import plotly.graph_objects as go
"""
data1 = pd.read_csv("multiOutput/summary/newData1114.tsv",
sep="\t",
dtype={'Pop':'str', 'Run':'int', 'Gen':'int', 'preAPref':'float', 'preIPref':'float', 'preOPref':'float', 'Matings':'float','Contacts':'float', 'MMContacts':'float', 'MalF':'float', 'FemF':'float', 'APref':'float', 'IPref':'float', 'OPref':'float','migrations':'float', 'A':'float', 'I':'float', 'O':'float', 'M':'float', 'F':'float', 'T':'float', 'matingSuccess':'float', 'misIdent':'float'})
data2 = pd.read_csv("multiOutput/summary/newData1115.tsv",
sep="\t",
dtype={'Pop':'str', 'Run':'int', 'Gen':'int', 'preAPref':'float', 'preIPref':'float', 'preOPref':'float', 'Matings':'float','Contacts':'float', 'MMContacts':'float', 'MalF':'float', 'FemF':'float', 'APref':'float', 'IPref':'float', 'OPref':'float','migrations':'float', 'A':'float', 'I':'float', 'O':'float', 'M':'float', 'F':'float', 'T':'float', 'matingSuccess':'float', 'misIdent':'float'})
data3 = pd.read_csv("multiOutput/summary/newData1118.tsv",
sep="\t",
dtype={'Pop':'str', 'Run':'int', 'Gen':'int', 'preAPref':'float', 'preIPref':'float', 'preOPref':'float', 'Matings':'float','Contacts':'float', 'MMContacts':'float', 'MalF':'float', 'FemF':'float', 'APref':'float', 'IPref':'float', 'OPref':'float','migrations':'float', 'A':'float', 'I':'float', 'O':'float', 'M':'float', 'F':'float', 'T':'float', 'matingSuccess':'float', 'misIdent':'float'})
data =pd.concat([data1, data2, data3])
"""
data1 = pd.read_csv("multiOutput/summary/newData1212.tsv",
sep="\t",
dtype={'Pop':'str', 'Run':'int', 'Gen':'int', 'preAPref':'float', 'preIPref':'float', 'preOPref':'float', 'Matings':'float','Contacts':'float', 'MMContacts':'float', 'MalF':'float', 'FemF':'float', 'APref':'float', 'IPref':'float', 'OPref':'float','migrations':'float', 'A':'float', 'I':'float', 'O':'float', 'M':'float', 'F':'float', 'T':'float', 'matingSuccess':'float', 'misIdent':'float'})
data2 = pd.read_csv("multiOutput/summary/newData1209.tsv",
sep="\t",
dtype={'Pop':'str', 'Run':'int', 'Gen':'int', 'preAPref':'float', 'preIPref':'float', 'preOPref':'float', 'Matings':'float','Contacts':'float', 'MMContacts':'float', 'MalF':'float', 'FemF':'float', 'APref':'float', 'IPref':'float', 'OPref':'float','migrations':'float', 'A':'float', 'I':'float', 'O':'float', 'M':'float', 'F':'float', 'T':'float', 'matingSuccess':'float', 'misIdent':'float'})

data =pd.concat([data1, data2])
data = data[pd.notnull(data["FFec"])]
data = data[data["Run"]%100 > 80]

data["tMorphs"] = data["A"]+data["I"]+data["O"]
data = data[data["tMorphs"]==1]

data.set_index(["Pop", "Run", "Gen"], inplace=True)
finalGen = max(data.index.get_level_values('Gen'))
finals = data.iloc[data.index.get_level_values('Gen')==finalGen]


morphs = np.array([finals["A"], finals["I"], finals["O"]])
print(morphs.shape)


fig = go.Figure(go.Scatterternary(a=data["A"], b=data["I"], c=data["O"], mode="markers"))
fig.show(renderer="browser")
