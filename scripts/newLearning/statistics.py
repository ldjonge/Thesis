import numpy as np
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from heatmap import *
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas_profiling
from scipy.optimize import curve_fit
from scipy import stats
import sys

def readFile(infile): # Read in the required datafile and prepare the data for analysis
    data = pd.read_csv(infile, sep="\t", dtype={'Pop':'str', 'Run':'int', 'Gen':'int', 'preAPref':'float', 'preIPref':'float', 'preOPref':'float', 'Matings':'float','Contacts':'float', 'MMContacts':'float', 'MalF':'float', 'FemF':'float', 'APref':'float', 'IPref':'float', 'OPref':'float','migrations':'float', 'A':'float', 'I':'float', 'O':'float', 'M':'float', 'F':'float', 'T':'float', 'matingSuccess':'float', 'misIdent':'float'}) # Read in required data file
    data.fillna(value=0, inplace=True)
    # data = data[data["Pop"] != "Total"]
    data.rename(index=str, columns={"T":"Total"}, inplace=True)
    data["pres"]=data["Total"].astype(bool).astype(int) # Allowing for easy removal of 'empty' populations from the analysis

    data["APres"] = (data["A"]>0).astype(int) # Measuring morphs presence/absence
    data["IPres"] = (data["I"]>0).astype(int)
    data["OPres"] = (data["O"]>0).astype(int)
    data["nMorphs"] = data["APres"]+data["OPres"]+data["IPres"] # Record number of morphs present


    return data


monomorphic = data[data["nMorphs"]==1] # For easy comparison between populations with different numbers of morphs present
dimorphic = data[data["nMorphs"]==2]
trimorphic = data[data["nMorphs"]==3]
