#!/usr/bin/env python3
import random

phenoDict = {"pp":"A", "pq":"A", "pr":"A", "qp":"A", "rp":"A", "qq":"I", "qr":"I", "rq":"I", "rr":"O"}

class Male:
    def __init__(self, pAll, mAll):
        self.genotype = pAll+mAll


class Female:
    def __init__(self,pAll,mAll):
        self.genotype = pAll+mAll
        self.phenotype = phenoDict[self.genotype]

def randomAllele(p,q,r):
    if random.random() <= p:
        allele = "p"
    elif random.random() <= q:
        allele = "q"
    else:
        allele = "r"
    return allele

#Population is always a list of Males and Females
def createPop(N, p, q, r):
    population = []
    pAll = randomAllele(p,q,r)
    mAll = randomAllele(p,q,r)
    for i in range(N):
        pAll = randomAllele(p,q,r)
        mAll = randomAllele(p,q,r)
        if i < N/2:
            population.append(Male(pAll,mAll))
        else:
            population.append(Female(pAll,mAll))
    return(population)


def reproduce(Male, Female, generation):
    if random.random() <= 0.5:
        pAll = Male.genotype[0]
    else:
        pAll = Male.genotype[1]
    if random.random() <= 0.5:
        mAll = Female.genotype[0]
    else:
        mAll = Female.genotype[1]
    if random.random() <= 0.5:
        generation.append(Male(pAll, mAll))
    else:
        generation.append(Female(pAll, mAll))

if __name__ == "__main__":
    print(createPop(10, 0.33, 0.33, 0.33))
