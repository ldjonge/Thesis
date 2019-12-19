# Discussion Points
## Migration
Much higher migration rates are required to form new populations than to introduce a new morph into an existing population. This in combination with NFDS maintaining diversity in any population that is big enough (which can still be relatively small), will lead to any new population forming quickly becoming at least dimorphic.

## Gynochrome frequency
All simulations have consistently shown that in populations where both gynochrome morphs are present, both morphs will stabilise at the same equilibrium frequency, because other than the dominance hierarchy they are identical in every way other than mating harassment. This means that the consistent difference in morph frequencies observed in real populations is evidence of a fitness difference between these morphs not directly related to learnt mate recognition. Whether this is a 'raw' fitness difference, a difference in optimal environment, or a difference in the detection by males obviously cannot be determined from these simulations but some difference must be present.


# Findings
## Differential maturation chance
It appears that under current conditions with no migration, all morphs are present at approximately equal frequency if gynochromes are 0.9 times as likely to reach maturity as androchromes and males. At 75% gynochrome maturation chance androchrome frequency appears to reach around the Skåne levels. If the difference gets bigger, populations tend to fix at a 100% androchrome frequency unless migration is introduced. Note that due to the way maturation chance is implemented, if androchromes and gynochromes have the same number of eggs, a 75% maturation chance of gynochromes leads to only a 57% androchrome frequency. 

## Migration
Migration should be modelled as only happening after the mating phase, or perhaps one 'migration session' before the mating phase and one after. Having migration occur more often increases the duration of the simulations too much. I'm thinking once before once after, this way genes of two populations can mix directly, and if a 'new' female enters a population before the mating phase she instantly has the rare morph advantage, whereas the same female entering fertilised after the mating phase means she does not yet have an advantage and as such the new morph is more likely to be eliminated by drift directly.

## Ternary Contour Plots
Still difficulty creating simulations which produce data covering the entire frequency space, working on it.
