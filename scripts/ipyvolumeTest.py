import numpy as np
import pandas as pd
import ipyvolume as ipv

data = pd.read_csv("multiOutput/summary/newData.tsv", sep="\t", dtype={'Pop':'str', 'Run':'int', 'Gen':'int', 'preAPref':'float', 'preIPref':'float', 'preOPref':'float', 'Matings':'float','Contacts':'float', 'MMContacts':'float', 'MalF':'float', 'FemF':'float', 'APref':'float', 'IPref':'float', 'OPref':'float','migrations':'float', 'A':'float', 'I':'float', 'O':'float', 'M':'float', 'F':'float', 'T':'float', 'matingSuccess':'float', 'misIdent':'float'})
#data.fillna(value=0, inplace=True)
ipv.figure(width=500, height=500, controls=True, controls_vr=True, debug=True)
ipv.xyzlabel("A","I","O")
ipv.scatter(data["A"], data["I"], data["O"], size=1, color="blue", marker="sphere")
ipv.style.set_style_dark()
ipv.save("morphDist.html")
