# Findings
A learning period before mating increases the probability of rare morphs going extinct by increasing morph preference to the point where rare morphs are at a significant risk of not mating.
While it is good at maintaining polymorphisms in larger populations or when closer to equilibrium frequencies, when present in small frequencies the effects of drift are enhanced

## Polymorphisms and Population Size
Having more morphs present in a population increases the average fitness of the females by reducing mating harassment, and as a result increases equilibrium population size. Monomorphic populations then have the lowest fitness, with mono-A populations having the lowest fitness due to their reduced fecundity not being counteracted by the effects of their mimicry.
Under current simulations trimorphic populations have an equilibrium size of 535, dimorphic 466 and monomorphic 317.
It should be noted that due to presumably both reduced baseline fitness and reduced male fitness populations with androchromes have lower equilibrium population sizes.
Dimorphic populations with A present have an equilibrium size of 454, those without A present 477.
Monomorphic populations are A:288, I:327 and O:336. The difference between I and O has a p-value of 0.012, all other p-values < 10^-10, correcting for the number of tests performed there is no significant difference between mono-I and mono-O populations. In subsequent simulations p-values for differences between A and H remained highly significant, the difference between I and O was insignificant.

## Introducing a new morph
*Androchromes tend to go extinct in trimorphic populations without interference, this influences these results*
So far no success in introducing a third morph of any kind
Effectivity of introducing a new morph depends on the morph based on two factors:
* Dominance Hierarchy: dominant alleles have a higher chance of staying 'relevant' in the population. This is logical, as introducing a homozygous recessive female into the population guarantees the next generation will not have this phenotype present, and thus relies on carriers to stay present for selection to potentially have an effect. Getting an introduction of Obsoleta into a population to lead to permanent presence appears to be near impossible.
* Baseline fecundity: While androchromes appear to be better at surviving the first few generations due to the allele being dominant, they are more likely to go extinct later as a result of their reduced baseline fecundity.
A learning phase massively reduces the chances of a new morph being successfully maintained in a population, probably as a result of the PFDS present for very rare morphs.
Introducing a single (unfertilised) female into a population has a relatively low chance of the morph being maintained. This can have multiple implications:
* Migration is unlikely to be effective unless a female is fertilised first
* Migration needs to happen in higher numbers to have permanent effects
* Migration happens repeatedly, and either a second introduction of the same morph happens before the morph has gone extinct, or due to enough repetition eventually one will be effective.

# Preliminary Conclusions
**Selection works on morph frequency, drift on allele frequency**
Maintaining a morph at a low frequency is easier if the allele for said morph is recessive, as this way the allele can be carried in other individuals. Equal morph frequency will occur at a higher allele frequency, so drift has a smaller chance of eliminating the morph and NFDS can keep it present at low frequencies.
Introducing a new morph however is much easier if the allele for the morph is dominant. If a single unfertilised female is introduced into a population and the allele for her morph is recessive, regardless of her own fitness there will be no expression of the allele in the next generation, even if many carriers are present. Therefore even if the entire next generation is made up of her offspring, only 1/4 of the females 2 generations later will express the new morph. In reality as there is no phenotypical effect of the allele in the carriers, complete random mating will occur, and the chances of 2 carriers mating are relatively low, with only 1/8 of the offspring (1/4 of the females) expressing the phenotype. Therefore extremely strong selective pressures are needed to maintain the new morph.

# To investigate
We know that the number of morphs affects population size, extinction risk is also based on the 'robustness' to random fluctuations. This robustness should be studied. Intuitively any option could make sense, populations with higher fecundity should be less prone to stochasticity, but it should be noted that random fluctuations in population size could randomly eliminate a morph, which could exacerbate the effects and lead to higher extinction risk in that way.
I can't find any publications about differential extinction risk in any species. While it makes intuitive sense to me, how to go about wording it and referencing the concept?
Perhaps link to specialisation, more specialised species are less resistant to changes, monomorphic could resemble specialisation in this case? Seems like a wonky argument at best.

# Figure out how to make ternary contour plots properly!!
