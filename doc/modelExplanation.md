# Individual Based Model
## Individuals
Each individual is a member of a Python Class
### Combined traits
Genotype is determined from parental alleles in $__init__$ function.  
Female phenotype is determined by the genotype using a dictionary to translate, male phenotype is set as male, but could be changed to investigate IASC at a later stage.
### Male
Fecundity determined in separate function, called after a full population is generated. Currently male fecundity is only depending on population size and carrying capacity, however it can start varying depending on environmental variables and possibly female phenotype frequency later
### Female
Fecundity determined in separate function, called after a full population is generated. Currently only dependent on population size, phenotype frequency and carrying capacity. Will later also depend on environmental variables.

## Populations
A first generation is randomly generated using the starting allele frequencies. To start there is an even 50/50 split between males and females.  
### Current model
For each consecutive generation each female has a chance to reproduce depending on her fecundity. If fecundity is above 1, she will produce at least one offspring, with the possibility of a second being (fecundity-1). A male is randomly sampled from the population, with sampling weighted based on male fecundity. Each offspring will have a 50\% chance to be either male or female, and the genotype consists of a random paternal and maternal allele, as in real reproduction.

### Potential model 1: Males' right to choose
As in reality males are not randomly sampled, rather males choose specific females based on learned preference, males may have to be the first parent to be chosen. The male may then first choose a colour depending on morph frequency, after which a specific individual female is chosen with that colour. Perhaps using something similar to \cite{Hardling2006}, letting the probability depend on morph frequencies.

### Potential model 2: Why stay in one place?
Based on the model with a single population, add in multiple populations that mostly evolve independently, but migration between the populations is possible at a relatively low frequency. Individuals that move may at first be sampled completely randomly, but can later be changed to be sampled based on sex and phenotype.

### Potential model 3: The new frontier
Expanding on the model with multiple populations, add in a more geographical model where populations are not all connected equally. Give 'border' populations a possibility for migration into a new area, where a new population is started. If the stars align a new population can be started, should be interesting to see how the phenotypes develop there.

### Model addition 1: Lay all the eggs
Rather than having one individual offspring produced by each mating, the more realistic option is to have many offspring generated. In this case each potential combination would make up 25% of the offspring from a single mating session (Punnett squares). A single copulation may for example produce 100 eggs in this scenario. Fecundity would no longer depend on carrying capacity in this case, rather a percentage of eggs hatches based on both the number of eggs each female lays and the carrying capacity. This model is more realistic and therefore more likely to account for drift in an accurate way.

### Model addition 2: It ain't easy being blue
Give andromorphs a lower baseline fecundity, accounting for the phenotypical changes other than colour that are associated with the A phenotype. This would of course be affected by geographical factors as well, probably leading to the northernmost populations still having andromorphs as the morph with the highest fecundity before sexual selection.

### Model addition 3: Choice matters
If males only search for the least common phenotype one would assume their mating success is reduced as a result. Chance of successful reproduction in males should be dependent on morph preference to 'encourage' good learning. This is probably accounted for already by having a preference for the most common morph, although there may need to be a chance of failing to find a mate if the total probability of mating is too low.

### Model addition 4: Do they ever learn?
While the evidence for learned preference is clear, it has not been excluded with certainty that genetics play no role in the mate preference of males. Having either the same locus or a different one play a role in mate preference could provide interesting clues about the potential interaction between genetics and plastic learning.

### Model addition 5: A little gay is okay
We could also add the possibility of a male with a preference for androchrome females accidentally choosing another male as mating partner, leading to a failed attempt that does not lead to reproduction. This would then result in the potential reduced success of males in a population with high androchrome frequencies. This accounts for the fact that apart from simply having multiple different morphologies, the androchrome females are in fact male mimics, which should have an additional advantage.

### Model addition 6: Leave me alone
It is known that females are negatively affected by mating harassment. This should obviously be accounted for in the model by reducing the fecundity. However as it stands the model does not contain a direct effect of mating on female fecundity. Presumably a good way to do this is reduce the female's fecundity every time she has been clasped.

### Model addition 7: Man knows man
High frequency of males may increase the preference for androchrome females.
