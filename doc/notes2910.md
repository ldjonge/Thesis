# Exclusions
Genetic basis for learning. Very interesting but currently too complicated to add, perhaps later but for now lets assume not. All males are phenotypically identical until they start learning.
This also means evolution of learning ability will not happen, all parameters including learning will be constant over the duration of any simulation, with the obvious exception of morph frequency.


# Things to implement, how?
Costs of learning/memory
Forgetting -- As preference is currently defined as (1+x), perhaps dividing x for each morph by 2? every sequence where no interaction occurs? Should probably increase the chance of successful mate finding if that is introduced as well, as otherwise preference will hardly change over many sequences.
Should female propensity to mate differ over time? Start low, increase until first mating, then reduce once mated? If so, should we assume a female knows that she is fertilised and will be more resistant permanently, or is there a possibility she 'forgets'?
Having virgin females be more promiscuous seems to solve a ton of problems beautifully. It seems to lead to potentially maintained trimorphism, with fecundity relatively equal between all morphs, suggesting this could be an equilibrium solution.
Importantly, this process intuitively fits with the idea that less common morphs seem to be more promiscuous, as these are less likely to mate early in their life cycle. As such they are on average virgins for a longer time, and thus will increase their willingness to mate.

# Things to implement, doable but not done yet
Simplify model, i.e. A vs H, I vs O only, this includes introducing a way that dimorphic populations are treated as such, and there is no constant preference for non-existent morphs. Adjusting this for the potential introduction of a new morph through migration will be interesting, but could almost certainly be done.
Monomorphic and Dimorphic populations now handle preference properly without being influenced by a non-existent morph, this way it will also make populations without androchromes easier to simulate due to the lack of potential misidentification. This probably means we need to recalculate morph frequencies after every possible migration moment. 

A vs H model may have to exclude learning altogether if we want to distinguish between mimicry and LMR, instead simply reducing harassment based on A frequency

Reduce probability of finding a mate based on presence of preferred morph, males won't actually choose between 200 females.
## Ideas:
Randomly sample 10? individuals from a population as potential mates, male chooses between these based on preference, if preferred morph isn't there than much reduced chance of mating.
Average preferences for each of the 10 individuals, this is chance of finding a partner, then weighted choice from the sample.
