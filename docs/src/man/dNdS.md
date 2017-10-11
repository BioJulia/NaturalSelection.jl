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
