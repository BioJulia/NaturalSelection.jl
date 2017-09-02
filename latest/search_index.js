var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#NaturalSelection.jl-1",
    "page": "Home",
    "title": "NaturalSelection.jl",
    "category": "section",
    "text": "(Image: Latest Release) (Image: NaturalSelection) (Image: License) (Image: BioJulia maintainer: bicycle1885) (Image: BioJulia maintainer: Ward9250)Development builds: (Image: Build Status) (Image: Build status) (Image: codecov)"
},

{
    "location": "index.html#Description-1",
    "page": "Home",
    "title": "Description",
    "category": "section",
    "text": "NaturalSelection.jl provides methods for detecting the presence, strength and effects of natural selection, in biological data."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Install NaturalSelection.jl from the Julia REPL:julia> Pkg.add(\"NaturalSelection\")If you are interested in the cutting edge of the development, please check out the master branch to try new features before release."
},

{
    "location": "man/dNdS.html#",
    "page": "dNdS",
    "title": "dNdS",
    "category": "page",
    "text": "CurrentModule = NaturalSelection.dNdS"
},

{
    "location": "man/dNdS.html#dN-/-dS-1",
    "page": "dNdS",
    "title": "dN / dS",
    "category": "section",
    "text": "The NaturalSelection.dNdS module provides several different methods of inferring the action of natural selection from coding sequences.Evolutionary pressures on proteins are often quantified by the ratio of substitution rates at non-synonymous and synonymous sites i.e. dN/dS.The dN/dS ratio was originally developed for application to distantly diverged sequences, the differences among which represent substitutions that have fixed along independent lineages.Nevertheless, the dN/dS measure is often applied to sequences sampled from a single population, the differences among which represent segregating polymorphisms. However, do be careful if this is what you are doing, as it has been demonstrated that dN/dS is not always suitable for such purposes (Sergey Kryazhimskiy & Joshua B. Plotkin, 2008)."
},

{
    "location": "man/dNdS.html#NaturalSelection.dNdS.dNdS_NG86",
    "page": "dNdS",
    "title": "NaturalSelection.dNdS.dNdS_NG86",
    "category": "Function",
    "text": "dNdS_NG86(x, y, k::Float64 = 1.0, code::GeneticCode)\n\nCompute dN and dS, using the Nei and Goborjei 1986 method.\n\nThis function requires two iterables x and y, which yield DNACodon or RNACodon type variables. These two types are defined in the BioSequences package.\n\n\n\ndNdS_NG86(x::BioSequence{A}, y::BioSequence{A}, k::Float64, code::GeneticCode) where {A <: NucAlphs}\n\nCompute dN and dS, using the Nei and Goborjei 1986 method.\n\nThis method adds conveinience when working with DNA or RNA sequences, by taking two sequences, and creating two vectors of aligned codons from them. These two iterables are then passed into the generic NG86 method.\n\n\n\n"
},

{
    "location": "man/dNdS.html#NaturalSelection.dNdS.S_N_NG86",
    "page": "dNdS",
    "title": "NaturalSelection.dNdS.S_N_NG86",
    "category": "Function",
    "text": "S_N_NG86(codon::C, k::Float64, code::GeneticCode) where {C <: CDN}\n\nEnumerate the number of synonymous (S) and non-synonymous (N) sites in a codon, using the method used by Nei and Goborjei (1986).\n\nReturns a tuple where S is the first element and N is the second (S, N).\n\nEach site in a codon may be both partially synonymous and non-synonymous.\n\n\n\n"
},

{
    "location": "man/dNdS.html#NaturalSelection.dNdS.DS_DN_NG86",
    "page": "dNdS",
    "title": "NaturalSelection.dNdS.DS_DN_NG86",
    "category": "Function",
    "text": "DS_DN_NG86(x::C, y::C, code::GeneticCode) where C <: CDN\n\nCompute the number of synonymous (DS) and non-synonymous (DN) mutations between two codons, using the all paths method used by the Nei and Goborjei (1986).\n\n\n\n"
},

{
    "location": "man/dNdS.html#The-NG86-method-1",
    "page": "dNdS",
    "title": "The NG86 method",
    "category": "section",
    "text": "dNdS_NG86\nS_N_NG86\nDS_DN_NG86"
},

{
    "location": "contributing.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing.html#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.If you have a question about contributing or using this package, you are encouraged to use the Bio category of the Julia discourse site.Detailed guidance for contributing to all BioJulia packages is provided at the BioJulia Contribution Documentation.Here we list specific details about contributing and maintainership pertaining specifically to the NaturalSelection.jl package."
},

{
    "location": "contributing.html#Named-maintainers-1",
    "page": "Contributing",
    "title": "Named maintainers",
    "category": "section",
    "text": "The named maintainers of this package is Ben Ward. It is their responsibility to make final choices about pull requests and issues, although because of our community structure, you will find other maintainers assisting them."
},

{
    "location": "contributing.html#Branching-model-1",
    "page": "Contributing",
    "title": "Branching model",
    "category": "section",
    "text": "The branching model used to develop and make releases of this package is the OneFlow model summarized in the BioJulia Contribution Documentation"
},

]}
