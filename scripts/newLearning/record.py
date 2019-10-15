from baseFunctions import *
import numpy as np
from pops import *

"""
Fecundity parameters should be recorded, not per individual but as a population average
These parameters are recorded per sex, as well as per female morph
Any morph or sex that is not present will have fecundity values set to NaN to make sure the data does not get corrupted.
"""
def recordFecStats(pop):
    if len(pop[0]) > 0:
        totalMFer = 0
        totalMMSucc = 0
        totalMSurv = 0
        for m in pop[0]:
            totalMFer += m.fertility
            totalMMSucc += m.mSucc
            totalMSurv += m.surv
        avgMFer = totalMFer/len(pop[0])
        avgMMSucc = totalMMSucc/len(pop[0])
        avgMSurv = totalMSurv/len(pop[0])
    else:
        avgMFer = np.nan
        avgMMSucc = np.nan
        avgMSurv = np.nan
    if len(pop[1]) > 0:
        totalFFec = 0
        totalFMSucc = 0
        totalFSurv = 0
        totalAFec = 0
        totalIFec = 0
        totalOFec = 0
        totalAMSucc = 0
        totalIMSucc = 0
        totalOMSucc =0
        totalASurv = 0
        totalISurv = 0
        totalOSurv = 0
        totalA = 0
        totalI = 0
        totalO = 0
        for f in pop[1]:
            totalFFec += f.fecundity
            totalFMSucc += f.mSucc
            totalFSurv += f.surv
            if f.phenotype == "A":
                totalAFec += f.fecundity
                totalAMSucc += f.mSucc
                totalASurv += f.surv
                totalA += 1
            elif f.phenotype == "I":
                totalIFec += f.fecundity
                totalIMSucc += f.mSucc
                totalISurv += f.surv
                totalI += 1
            elif f.phenotype == "O":
                totalOFec += f.fecundity
                totalOMSucc += f.mSucc
                totalOSurv += f.surv
                totalO += 1
        avgFFec = totalFFec/len(pop[1])
        avgFMSucc = totalFMSucc/len(pop[1])
        avgFSurv = totalFSurv/len(pop[1])
        if totalA > 0:
            avgAFec = totalAFec/totalA
            avgAMSucc = totalAMSucc/totalA
            avgASurv = totalASurv/totalA
        else:
            avgAFec = np.nan
            avgAMSucc = np.nan
            avgASurv = np.nan
        if totalI > 0:
            avgIFec = totalIFec/totalI
            avgIMSucc = totalIMSucc/totalI
            avgISurv = totalISurv/totalI
        else:
            avgIFec = np.nan
            avgIMSucc = np.nan
            avgISurv = np.nan
        if totalO > 0:
            avgOFec = totalOFec/totalO
            avgOMSucc = totalOMSucc/totalO
            avgOSurv = totalOSurv/totalO
        else:
            avgOFec = np.nan
            avgOMSucc = np.nan
            avgOSurv = np.nan
    else:
        avgFFec = np.nan
        avgFMSucc = np.nan
        avgFSurv = np.nan
        avgAFec = np.nan
        avgAMSucc = np.nan
        avgASurv = np.nan
        avgIFec = np.nan
        avgIMSucc = np.nan
        avgISurv = np.nan
        avgOFec = np.nan
        avgOMSucc = np.nan
        avgOSurv = np.nan
    return [avgMFer, avgMMSucc, avgMSurv, avgFFec, avgFMSucc, avgFSurv, avgAFec, avgAMSucc, avgASurv, avgIFec, avgIMSucc, avgISurv, avgOFec, avgOMSucc, avgOSurv]

"""
Preference should be recorded, again as a population average.
Again in populations with no males the preferences are set to NaN.
"""
def recordPref(pop):
    totalAPref = 0
    totalIPref =0
    totalOPref = 0
    for m in pop[0]:
        prefSm = sum(m.prefs.values())
        totalAPref += m.prefs["A"] #(m.prefs["A"]/prefSm*0.7+0.1)
        totalIPref += m.prefs["I"] #(m.prefs["I"]/prefSm*0.7+0.1)
        totalOPref += m.prefs["O"] #(m.prefs["O"]/prefSm*0.7+0.1)
    prefSum = totalAPref + totalIPref + totalOPref
    if prefSum > 0:
        avgAPref = totalAPref/len(pop[0])
        avgIPref = totalIPref/len(pop[0])
        avgOPref = totalOPref/len(pop[0])
        return [avgAPref, avgIPref, avgOPref]
    else:
        return [np.nan,np.nan,np.nan]

def preRecord(freqTable, pops, gen):
    for pop in pops:
        id = pops.index(pop)+1
        prefs = recordPref(pop)
        phenFreq = calcPhenoFreq(pop)
        freqTable.append([id, gen, "preAPref", prefs[0]])
        freqTable.append([id, gen, "preIPref", prefs[1]])
        freqTable.append([id, gen, "preOPref", prefs[2]])
        freqTable.append([id, gen, "A", phenFreq["A"]])
        freqTable.append([id, gen, "I", phenFreq["I"]])
        freqTable.append([id, gen, "O", phenFreq["O"]])
        freqTable.append([id, gen, "M", len(pop[0])])
        freqTable.append([id, gen, "F", len(pop[1])])
        freqTable.append([id, gen, "T", len(pop[0])+len(pop[1])])

def postRecord(freqTable, pops, matings, contacts, deaths, gen):
    for pop in pops:
        id = pops.index(pop)
        avgFecs = recordFecStats(pop)
        prefs = recordPref(pop)
        #print(id, prefs)
        freqTable.append([id+1, gen, "MFer", avgFecs[0]])
        freqTable.append([id+1, gen, "MMSucc", avgFecs[1]])
        freqTable.append([id+1, gen, "MSurv", avgFecs[2]])
        freqTable.append([id+1, gen, "FFec", avgFecs[3]])
        freqTable.append([id+1, gen, 'FMSucc', avgFecs[4]])
        freqTable.append([id+1, gen, 'FSurv', avgFecs[5]])
        freqTable.append([id+1, gen, 'AFec', avgFecs[6]])
        freqTable.append([id+1, gen, 'AMSucc', avgFecs[7]])
        freqTable.append([id+1, gen, 'ASurv', avgFecs[8]])
        freqTable.append([id+1, gen, 'IFec', avgFecs[9]])
        freqTable.append([id+1, gen, 'IMSucc', avgFecs[10]])
        freqTable.append([id+1, gen, 'ISurv', avgFecs[11]])
        freqTable.append([id+1, gen, 'OFec', avgFecs[12]])
        freqTable.append([id+1, gen, 'OMSucc', avgFecs[13]])
        freqTable.append([id+1, gen, 'OSurv', avgFecs[14]])
        freqTable.append([id+1, gen, 'Matings', matings[id]])
        freqTable.append([id+1, gen, 'Contacts', contacts[id]])
        #freqTable.append([id+1, gen, 'MMContacts', MMcontacts[id]])
        #freqTable.append([id+1, gen, 'Migrations', migrations[id]])
        #freqTable.append([id+1, gen, 'Fertilised Migrations', fertMigrations[id]])
        freqTable.append([id+1, gen, 'Deaths', deaths[id]])
        freqTable.append([id+1, gen, "APref", prefs[0]])
        freqTable.append([id+1, gen, "IPref", prefs[1]])
        freqTable.append([id+1, gen, "OPref", prefs[2]])
