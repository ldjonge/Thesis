import random
from pops import *
from numpy.random import choice


"""
Males will search for mates based on their preference.
The number of interactions per population/generation should be recorded.
"""
def matingSearch(pop, params, popDict):
    matings = 0
    contacts = 0
    #MMcontacts = 0
    deaths = 0
    totalPop = pop[0] + pop[1]
    N = len(totalPop)
    if len(pop[1])!=0:
        for male in pop[0]:
            if male.taken == 0: #Copulating males cannot search for a new mate
                # the chance of approaching a specific individual is determined for each potential pair, based on the male's preference.
                hitChance = []
                male.calcmSucc()
                """
                for ind in pop[0]:
                    # In case there should be the potential for a male to approach another male
                    if ind.taken != 0:
                        hitChance.append(0)
                    else:
                        hitChance.append(max(params["maleRec"]*male.prefs["A"], 0))
                """
                for fem in pop[1]:
                    if fem.taken != 0:
                        hitChance.append(0)
                    elif fem.phenotype == "A":
                        hitChance.append(max(male.prefs["A"]*fem.vis, 0))
                    elif fem.phenotype == "I":
                        hitChance.append(max(male.prefs["I"]*fem.vis, 0))
                    elif fem.phenotype == "O":
                        hitChance.append(max(male.prefs["O"]*fem.vis, 0))
                totalW = sum(hitChance)
                # Probabilities of selecting a specific individual are calculated from the weights
                if totalW != 0:
                    for i in range(len(hitChance)):
                        hitChance[i] = hitChance[i]/totalW
                    if random.random() < male.mSucc*N/popDict["K"]:    # Male mating success is used as probability of finding a mate
                        mate = choice(pop[1], p=hitChance)
                    else:
                        mate = None
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
                            mate.mSucc *= params["mateFFertEff"]
                            male.fertility *= params["mateMFertEff"]
                            mate.surv *= params["mateSurvEff"]
                            male.surv *= params["mateSurvEff"]
                            if mate.phenotype == "A":
                                male.prefs["A"] *=params["matePrefEff"]
                            elif mate.phenotype == "I":
                                male.prefs["I"] *=params["matePrefEff"]
                            elif mate.phenotype == "O":
                                male.prefs["O"] *=params["matePrefEff"]
                        else:
                            mate.fecundity *= params["failFecEff"]
                            # mate.mSucc *=params["failFertEff"]
                            mate.surv *= params["failSurvEff"]
                            if mate.phenotype == "A":
                                male.prefs["A"] *=params["failPrefEff"]
                            elif mate.phenotype == "I":
                                male.prefs["I"] *=params["failPrefEff"]
                            elif mate.phenotype == "O":
                                male.prefs["O"] *=params["failPrefEff"]
                    elif type(mate) == Male:
                        MMcontacts += 1
                        mate.fertility *= params["failFertEff"]
                        mate.surv *= params["failSurvEff"]
                        male.prefs["A"] *=params["failPrefEff"]
            else:
                male.taken -= 1
        #Individuals that die will be removed from the population
        pop[0] = [i for i in pop[0] if i.surv > random.random()]
        pop[1] = [i for i in pop[1] if i.surv > random.random()]
        deaths += (N - len(pop[0]) - len(pop[1])) # Number of deaths is recorded
    else:
        pass
    return [matings, contacts, deaths]

def migrate(pops, migrationMatrix):
    newPops = [([],[]) for i in pops]
    migrationList = []
    for pop in pops:
        migrations = 0
        fertMigrations = 0
        id = pops.index(pop)
        mRates = migrationMatrix.iloc[id]
        mChance = mRates.cumsum()
        for male in pop[0]:
            if male.migrate == 0:	#Individuals that are still migrating cannot migrate again
                migOut = random.random()
                for i in range(len(mChance)):
                    if migOut <= mChance[i]:
                        newPop = i
                        break
                if newPop != id:
                    migrations += 1
                    male.migrate += random.randint(1,5)
                    """
                    Depending on distance and speed migration may take a variable amount of time.
                    This could also be related to the migration probability depending on how that is decided on
                    """
                newPops[newPop][0].append(male)	#Individuals are added to the population for the next cycle
            else:
                male.migrate -= 1	# Migrating individuals get one step closer to their destination
                newPops[id][0].append(male)	# But for the purpose of recording are 'staying' in the population they are going towards
        for female in pop[1]:
            if female.migrate == 0:
                migChance = random.random()
                for i in range(len(mChance)):
                    if migOut <= mChance[i]:
                        newPop = i
                        break
                if newPop != id:
                    migrations += 1
                    female.migrate += random.randint(1,5)
                    if len(female.mates) > 0:
                        fertMigrations += 1
                newPops[newPop][0].append(female)
            else:
                female.migrate -= 1
                newPops[id][0].append(female)
        migrationList.append([migrations, fertMigrations])
    pops = newPops
    return migrationList
