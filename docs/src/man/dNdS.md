```@meta
CurrentModule = NaturalSelection
```

# dN / dS

NaturalSelection.jl provides several different methods of inferring
the action of natural selection from coding sequences.

Evolutionary pressures on proteins are often quantified by the ratio of
substitution rates at non-synonymous and synonymous sites i.e. dN/dS.

The dN/dS ratio was originally developed for application to distantly diverged
sequences, the differences among which represent substitutions that have fixed
along independent lineages.

Nevertheless, the dN/dS measure is often applied to sequences sampled from a
single population, the differences among which represent segregating
polymorphisms. However, do be careful if this is what you are doing, as it has
been demonstrated that dN/dS is not always suitable for such purposes
[(Sergey Kryazhimskiy & Joshua B. Plotkin, 2008)](https://doi.org/10.1371/journal.pgen.1000304).


## The NG86 method

```@docs
dNdS_NG86
S_N_NG86
DS_DN_NG86
```


# The MacDonald Kreitman Test

This test detects and measure the amount of adaptive evolution within a species
by determining whether adaptive evolution has occurred, and the proportion of
substitutions that resulted from positive selection.

To do this, the [McDonald–Kreitman](https://doi.org/10.1038%2F351652a0)
test compares the amount of variation within a
species (polymorphism) to the divergence between species (substitutions) at two
types of sites, neutral and nonneutral.

A substitution refers to a nucleotide that is fixed within one species, but a
different nucleotide is fixed within a second species at the same base pair of
homologous DNA sequences.

The two types of sites can be either synonymous or nonsynonymous within a
protein-coding region.

The null hypothesis of the McDonald–Kreitman test is that the ratio of
nonsynonymous to synonymous variation within a species is going to equal the
ratio of nonsynonymous to synonymous variation between species (i.e. Dn/Ds = Pn/Ps).

When positive or negative selection (natural selection) influences nonsynonymous
variation, the ratios will no longer equal. The ratio of nonsynonymous to
synonymous variation between species is going to be lower than the ratio of
nonsynonymous to synonymous variation within species (i.e. Dn/Ds < Pn/Ps) when
negative selection is at work, and deleterious mutations strongly affect
polymorphism. The ratio of nonsynonymous to synonymous variation within species
is lower than the ratio of nonsynonymous to synonymous variation between species
(i.e. Dn/Ds > Pn/Ps) when we observe positive selection.

Under neutrality the expectation is (Pn / Ps) == (Dn / Ds).

The statistic [α](https://doi.org/10.1038%2F4151022a) represents the proportion of
substitutions driven by positive selection. α can be equal to any number between
-Inf and 1.
Negative values of alpha are produced by sampling error, assumption violations,
or the segregation of [slightly deleterious](https://doi.org/10.1093/molbev/msn005)
amino acid mutations.
The null hypothesis here is that α = 0.

The neutrality index (NI) quantifies the direction and degree of departure from
neutrality (where Pn/Ps and Dn/Ds ratios equal).
A neutrality index greater than 1 is indicative of negative selection, a
neutrality index lower than 1 indicates positive selection is at work in the
population.

This test is provided by the function `mkt`.

```@docs
mkt
```

## References

For more about the McDonald Kreitman Test, see the following references:

https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test

[McDonald, J. H. Kreitman (1991)](https://doi.org/10.1038%2F351652a0)

[Charlesworth, J. Eyre-Walker (2008)](https://doi.org/10.1093%2Fmolbev%2Fmsn005)

[Eyre-Walker, A (2006)](http://www.zoology.wisc.edu/courses/611/Part2/Readings/3eyre-walker_tree2006.pdf)

[Stoletzki, N. Eyre-Walker (2010)](https://doi.org/10.1093%2Fmolbev%2Fmsq249)
