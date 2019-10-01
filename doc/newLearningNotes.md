# New Learning Model
The current learning model is based on the assumption that males may or may not have a preference upon emerging, which if present would be genetically determined. The possibility of learning before the first mating attempts should still be discussed, but either way interactions with females have an effect on mate preference. The preference for a morph will increase upon successful mating with said morph, and decrease after a failed mating attempt.

Mating attempts have a negative effect on female fecundity. It might be good to model in such a way that while the effect of mating attempts are the same whether successful or not, successful mating 'removes' the female from the mating population, and as such reduces the probability of future harassment.

All individuals in a generation are present at the same time. While not entirely realistic, it should hopefully not make too much of a difference.

## Males
Currently do not forget, possibly bring preference closer to even after every few searches?
### Parameters
* Innate preference (possibly based on genotype)
* Learning speed
* Survival?
* Fertility?

## Females
Can mate and be fertilised multiple times, but only lay eggs once at the end of a 'season'
In case of migration, this means a female could migrate both before and after fertilisation, but always before laying eggs.
'Father' of the offspring is determined based on the order of mating, the last male to mate will father 80% of the offspring, the rest will be evenly divided between all previous males. Of course if a female has only mated once, that one male will be the father of all the offspring.
### Parameters
* Mating resistance, should make sure that succChance*succEff > failChance*failEff
* Fecundity
* Survival?
*
