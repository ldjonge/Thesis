import random
import numpy as np
import pandas as pd
import sys
#import matplotlib.pyplot as plt
#import seaborn as sns; sns.set()
from readFiles import *
from baseFunctions import *
from record import *
from pops import *
from sims import *
import sampleMating

def runSim():
    paramDict = readParams()
    params = readPopInfo()
    pops = []
    nPop = len(params)
    for i in range(nPop):
        pops.append(startingPop(params[i], paramDict))
    freqTable = [["Pop", "Gen", "Pheno", "Value"]]
    for gen in range(1, int(paramDict["nGen"])+1): # First generation should be 1, not 0
        totalLen = 0
        matings = [0 for i in range(nPop)]
        contacts = [0 for i in range(nPop)]
        deaths = [0 for i in range(nPop)]
        popSizes = []
        freqs = preRecord(freqTable, pops, gen)
        for pop in pops:
            popSizes.append(len(pop[0])+len(pop[1]))
        for i in range(int(paramDict["genLength"])):
            for pop in pops:
                id = pops.index(pop)
                pres = [phen for phen in freqs[id].keys() if freqs[id][phen]!=0]
                results = sampleMating.matingSearch(pop, paramDict, params[id], pres)
                matings[id] += results[0]
                contacts[id] += results[1]
                deaths[id] += results[2]
                for fem in pop[1]:
                    if fem.taken != 0:
                        fem.taken -= 1
                for male in pop[0]:
                    male.calcmSucc()
        newPops = []
        totalLen = 0
        postRecord(freqTable, pops, matings, contacts, deaths, gen, freqs)
        for pop in pops:
            id = pops.index(pop)
            newPop = []
            #unMated = {"A":0, "I":0, "O":0}
            for f in pop[1]:
                f.eggLay(newPop, params[id])
            #    if len(f.mates) ==0:
            #        unMated[f.phenotype] += 1
            #print(len(pop[1]))
            #print(unMated)
            size = newPopSize(newPop, params[id]["K"])
            if gen%24==0:
                size = size/4
            newPop = popControl(newPop, size)
            popSize = len(newPop[0]) + len(newPop[1])
            if popSize > 0 and popSizes[id] ==0:
                print("New Population formed at {} in generation {}".format(id+1, gen))
            elif popSize == 0 and popSizes[id] > 0:
                print("Population {} extinct in generation {}".format(id+1, gen))
            totalLen += popSize
            #if gen == 50:
            #    for i in range(5):
            #        newPop[1].append(Female("r", "r"))
            phenFreq = calcPhenoFreq(newPop, male=True)
            for ind in newPop[1]:
                ind.calcFec(popInfo=params[id])
                ind.calcVis(phenFreq)
            for ind in newPop[0]:
                ind.calcFec(popInfo=params[id])
                #ind.learning(newPop[1], paramDict)
            newPops.append(newPop)
        if gen%10 == 0:
            print("Generation {} complete, population size {}".format(str(gen), str(totalLen)))
        pops = newPops
    return freqTable


if __name__ == "__main__":
    freqTable = pd.DataFrame(runSim())
    if len(sys.argv) > 1:
        freqTable.to_csv(sys.argv[1], header=None, index=None)
    else:
        pass
        # print(freqTable)
