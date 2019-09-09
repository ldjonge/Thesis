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
        popSizes = []
        preRecord(freqTable, pops, gen)
        for pop in pops:
            popSizes.append(len(pop[0])+len(pop[1]))
        for i in range(int(paramDict["genLength"])):
            for pop in pops:
                id = pops.index(pop)
                results = matingSearch(pop, paramDict, params[id])
                matings[id] += results[0]
                contacts[id] += results[1]
                for fem in pop[1]:
                    if fem.taken != 0:
                        fem.taken -= 1
        newPops = []
        totalLen = 0
        postRecord(freqTable, pops, matings, contacts, gen)
        for pop in pops:
            id = pops.index(pop)
            newPop = []
            for f in pop[1]:
                f.eggLay(newPop, paramDict["nEggs"], params[id])
            size = newPopSize(paramDict["nEggs"], newPop, params[id]["K"])
            newPop = popControl(newPop, size)
            popSize = len(newPop[0]) + len(newPop[1])
            if popSize > 0 and popSizes[id] ==0:
                print("New Population formed at {} in generation {}".format(id+1, gen))
            elif popSize == 0 and popSizes[id] > 0:
                print("Population {} extinct in generation {}".format(id+1, gen))
            totalLen += popSize
            for ind in newPop[1]:
                ind.calcFec(popInfo=params[id])
            for ind in newPop[0]:
                ind.calcFec(popInfo=params[id])
                ind.learning(newPop[1], paramDict)
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
        print(freqTable)
