# 'Logistics'
Finishing in Jan/Feb, how does this work given Costa Rica?
PhD possibility automatically follows as a discussion point

Half time should probably be sent in soon now, check what is needed and if this needs to be discussed at all

# S-shaped curves
## a) Steepness
In testing  
Steepness will probably be similar for all morphs, as there should not be any reason to assume a different feedback from different morphs. Successful mating is successful mating and failed mating is failed mating, unless there is some very clear evidence that some morph may give a bigger positive or negative feedback. Negative feedback could be interesting to study given behavioural differences, but may be too specific for not enough result??

## b) Inflection Point
In testing  
This is where morph differences will occur, as mimicry and crypticity will play a role here. As fitness is currently directly (negatively) correlated with male preference, a higher inflection point does indeed lead to a higher 'optimal' frequency.

## c) Intercept
Seems like a very realistic option, am trying it out right now. Should be less important if there is no learning pre-mating, however is still good to include, especially if we are to compare systems with and without learning pre-mating.

## d) Asymptote
Coincides with intercept, as the preference values as they are now always add up to 1 for the three morphs. Therefore the asymptote will be 1-(other morph intercepts)

# Other things to discuss
## Potential different way of modelling male fitness
While the Buridan's ass principle is an interesting option, it may also be the case that males get better at finding mates by being successful before. They may start at a certain medium success rate, as they are inexperienced and may not know very well what to do. If they do find a female and successfully mate, this should make them more experienced and 'confident' for finding another mate. If on the other hand they are unsuccessful in mating, or especially if they find a male, their mating success may go down as they are less 'confident'. For the most part this should lead to similar results as using preference as a proxy, but is perhaps more intuitive. In addition it also makes the idea that monomorphic populations have lower fitness potentially more real rather than created on purpose.

## Many of the current results may still be due to small population sizes
Mimicry in itself can currently counteract the effect of reduced fecundity in androchromes, but it should be looked at if androchromes go extinct in large populations when mimicry is not accounted for, NFDS by itself may do enough if populations are large enough. Mimicry may however be able to increase the equilibrium frequency disregarding

## Learning pre-mating
Currently it appears as if a learning phase before mating increases the effects of PFDS rather than NFDS, by making rare morphs almost impossible to be mated. However testing with intercepts above 0 may change this. Again should be noted that due to small population sizes a 3% frequency may in fact mean 3 individuals, drift will have a big effect.

## Population size/extinction risk
Extinction risk ~=~ Population size + Random variation in population size
In principle population size should be a clear indicator of extinction risk, smaller populations are more vulnerable to random shifts in population size. Currently not too much randomness in the way population sizes are determined for the next generation because in combination with small population sizes this increased the effects of drift even more by effectively introducing bottlenecks, however in later simulations randomness will be introduced and over multiple runs of the same simulations extinction risk can be evaluated.  
While population size would be the clearest indicator of extinction risk, it may still be good to indeed see if there is a difference in extinction risk introduced by varying the steepness of the S-curve due to potentially decreased/increased randomness in the population size.

Early indications that a steeper curve does decrease population size as the more common morph is targeted more. Stronger preference also increases the rates of males finding females, which additionally reduces the female fitness even further.

## Evolving the 'curves'
In nature they almost certainly can evolve and probably have evolved, however doing this in the simulations may not be worth focusing on. On the one hand it probably makes things more complicated, additionally the speed of evolution should be considered. Currently running most simulations with a duration of 100 generations, this is almost certainly not enough for any 'relevant' change to occur. For relevant change simulation length should probably be increased to 10,000 generations at least, this will take ages to run and may not be worth it if simulations with 1000 generations show mostly similar results to sims with 100 generations. One can imagine that no S-shaped curve was present in a monomorphic state, as there were no morph frequencies yet. The effects of a third morph on the curves of the previous 2 are interesting, although I could personally easily imagine that preference is based only on the presence of the morph itself rather than any effects of the other morph frequencies. Nonetheless if the tail-spot plays a big role in identifying females gynochromes could certainly be affected by each other's frequencies as well. Again interesting to discuss but possibly not worth modeling without clear evidence from the field.

# Analysing final data
Running multiple simulations and doing statistical analysis seems like the most reasonable solution for a highly stochastic model, no real way to approach it more analytically.
## What should be recorded?
Population size and morph frequencies are obvious.
Morph preference seems like an important analytical tool for the time being, but is probably not that important in the final analysis, though will probably show the S-curve which might be good
Male mating success could also be interesting to see if this is indeed correlated with highest preference, assuming we do not count it as a direct consequence.
Female fecundity should be important enough to consider, possibly other fitness components as well. 
Matings/contacts/deaths are probably not relevant for final results, migration should probably go in this as well. While all these numbers could be interesting in theory, they will generally be a direct result of input parameters.
