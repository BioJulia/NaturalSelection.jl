```@meta
CurrentModule = NaturalSelection
```

# Tajima's D

Tajima's D is a population genetic test statistic created by and named after the Japanese researcher Fumio Tajima.

Tajima's D is computed as the difference between two measures of genetic diversity:
The mean number of pairwise differences and the number of segregating sites, each scaled
so that they are expected to be the same in a neutrally evolving population of constant size.

The purpose of the statistic is to distinguish between a DNA sequence evolving
randomly ("neutrally") and one evolving under a non-random process.
The non-random process might be directional or balancing selection,
demographic expansion or contraction, genetic hitchhiking, or even introgression.

```@docs
tajimad(::AbstractFloat, ::Integer, ::Integer)
tajimad(::Any)
```
