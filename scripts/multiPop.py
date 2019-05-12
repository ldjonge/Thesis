import createPop
import pandas
from mating import *
from createPop import *
from standardFunctions import *

def migrate(pops, female):
    if female in pops[0][1]:
        pops[0][1].remove(female)
        pops[1][1].append(female)
    elif female in pops[1][1]:
        pops[1][1].remove(female)
        pops[0][1].append(female)

def main(length):
    paramDict = readParams()
    pop1 = startingPop(int(paramDict["N"]), paramDict["p"], paramDict["q"], paramDict["r"])
    pop2 = startingPop(int(paramDict["N"]), paramDict["p"], paramDict["q"], paramDict["r"])
    pops = [pop1, pop2]
    nEggs = paramDict["nEggs"]
    if length=="short":
        freqTable = [["Pop", "Gen", "A", "I", "O", "Males", "Females", "Total", "MaleFec", "FemFec", "APref", "IPref", "OPref", "Matings", "Contacts", "MMContacts"]]
    elif length=="long":
        freqTable = [["Pop", "Gen", "Pheno", "Value"]]
    for gen in range(1, int(paramDict["nGen"])+1):
        for pop in pops:
            id = pops.index(pop)+1
            nextGen = []
            matings = 0
            contacts = 0
            MMcontacts = 0
            if length == "long":
                prefs = recordPref(pop)
                freqTable.append([id, gen, "preAPref", prefs[0]])
                freqTable.append([id, gen, "preIPref", prefs[1]])
                freqTable.append([id, gen, "preOPref", prefs[2]])
            for i in range(30):
                results = matingSearch(pop, nEggs, paramDict["successRate"], nextGen, paramDict["K"])
                matings += results[0]
                contacts += results[1]
                MMcontacts += results[2]
                for fem in pop[1]:
                    if fem.taken != 0:
                        fem.taken -= 1
        for fem in pop1[1]:
            if random.random() < float(paramDict["m"]):
                migrate(pops, fem)
        for fem in pop2[1]:
            if random.random() < float(paramDict["m"]):
                migrate(pops, fem)
        for pop in pops:
            id = pops.index(pop)+1
            for fem in pop[1]:
                fem.eggLay(nextGen, nEggs)
            size = newPopSize(nEggs, nextGen, int(paramDict["K"]))
            avgFecs = recordFec(pop)
            prefs = recordPref(pop)

            pop = popControl(nextGen, size)
            totalPop = pop[0]+pop[1]
            phenFreq = calcPhenoFreq(pop)
            for ind in pop[0]:
                ind.calcFec(phenFreq)
                ind.complexLearning(totalPop)
            for ind in pop[1]:
                ind.calcFec(phenFreq)
            #print([str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1]))])
            if length == "short":
                freqTable.append([str(id), str(gen), str(phenFreq["A"]), str(phenFreq["I"]), str(phenFreq["O"]), str(len(pop[0])), str(len(pop[1])), str(len(pop[0])+len(pop[1])), str(avgFecs[0]), str(avgFecs[1]), str(prefs[0]), str(prefs[1]), str(prefs[2]), str(matings), str(contacts), str(MMcontacts)])
            elif length == "long":
                freqTable.append([id, gen, "A", phenFreq["A"]])
                freqTable.append([id, gen, "I", phenFreq["I"]])
                freqTable.append([id, gen, "O", phenFreq["O"]])
                freqTable.append([id, gen, "M", len(pop[0])])
                freqTable.append([id, gen, "F", len(pop[1])])
                freqTable.append([id, gen, "T", len(pop[0])+len(pop[1])])
                freqTable.append([id, gen, "MalF", avgFecs[0]])
                freqTable.append([id, gen,"FemF", avgFecs[1]])
                freqTable.append([id, gen, "APref", prefs[0]])
                freqTable.append([id, gen, "IPref", prefs[1]])
                freqTable.append([id, gen, "OPref", prefs[2]])
                freqTable.append([id, gen, "Matings", matings])
                freqTable.append([id, gen, "Contacts", contacts])
                freqTable.append([id, gen, "MMContacts", MMcontacts])

            else:
                print("failure")
        if (gen)%10 == 0:
            print("Generation {} complete, population size {}".format(str(gen), str(len(pops[0][0])+len(pops[0][1])+len(pops[1][1])+len(pops[1][0]))))

    return freqTable

if __name__ == "__main__":
    freqTable = main("short")
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'w') as outfile:
            for line in freqTable:
                print("\t".join(line), file=outfile)
    else:
        for line in freqTable:
            print("\t".join(line))
