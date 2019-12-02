import random
from pops import *
from numpy.random import choice
from baseFunctions import *

"""
Males will search for mates based on their preference.
The number of interactions per population/generation should be recorded.
"""
def matingSearch(pop, params, popDict, pres):
    attempts = 0
    matings = 0
    contacts = 0
    MMcontacts = 0
    deaths = 0
    totalPop = pop[0] + pop[1]
    N = len(totalPop)
    nonMating = []
    for male in pop[0]:
        for i in male.prefs.keys():
            if i not in pres:
                male.prefs[i] = 0
    for ind in totalPop:
        if ind.taken == 0:
            nonMating.append(ind)
    nPhens = len(pres)
    sampleSize = randomRound(10*len(nonMating)/popDict["K"])
    if len(pop[1])!=0:
        for male in pop[0]:
            if male.taken == 0:
                prefSum = sum(male.prefs.values())+nPhens
                if prefSum !=0:
                    attempts += 1
                    sample = random.sample(nonMating, sampleSize)
                    while len(sample) < 10:
                        sample.append(None)
                    hitChance = []
                    for ind in sample:
                        if type(ind)==Male:
                            hitChance.append(max(params["maleRec"]*((1+male.prefs["A"])/prefSum*(1-0.1*nPhens)+0.1*("A" in pres)), 0))
                        elif type(ind)==Female:
                            hitChance.append(max((1+male.prefs[ind.phenotype])/prefSum*(1-0.1*nPhens)+0.1*(ind.phenotype in pres),0))
                        else:
                            hitChance.append(0)
                    totalW = sum(hitChance)
                    mSucc = totalW/len(hitChance)
                    if totalW != 0:
                        for i in range(len(hitChance)):
                            hitChance[i] = hitChance[i]/totalW
                        if random.random() < mSucc:
                            mate = choice(sample, p=hitChance)
                        else:
                            mate=None
                        if type(mate)==Female:
                            contacts += 1
                            """
                            Female mating success incorporates both the chance of contact leading to mating and the chance of mating leading to fertilisation.
                            While biologically very different processes, the end result is very similar. Without fertilisation the female is harassed,
                            but no offspring is produced.
                            """
                            if random.random() <= mate.mSucc:
                                matings += 1
                                mate.mate(male, params)
                                mate.mSucc = popDict["{}MSucc".format(mate.phenotype)]
                                male.fertility *= params["mateMFertEff"]
                                mate.surv *= params["mateSurvEff"]
                                male.surv *= params["mateSurvEff"]
                                male.prefs[mate.phenotype] *= params["matePrefEff"]
                            else:
                                mate.fecundity *= params["failFecEff"]
                                mate.surv *= params["failSurvEff"]
                                male.prefs[mate.phenotype] *= params["failPrefEff"]
                        elif type(mate) == Male:
                            MMcontacts += 1
                            mate.fertility *= params["failFertEff"]
                            mate.surv *= params["failSurvEff"]
                            male.prefs["A"] *=params["failPrefEff"]
            else:
                male.taken -= 1
        #Individuals that die will be removed from the population
        #pop[0] = [i for i in pop[0] if i.surv > random.random()]
        #pop[1] = [i for i in pop[1] if i.surv > random.random()]
        for fem in pop[1]:
            if fem.taken==0:
                fem.mSucc *= 1.1
        deaths += (N - len(pop[0]) - len(pop[1])) # Number of deaths is recorded
    else:
        pass
    return [matings, contacts, deaths]
