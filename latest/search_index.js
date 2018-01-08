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
    "text": "Latest release:(Image: Latest Release) (Image: NaturalSelection) (Image: License) (Image: ) (Image: BioJulia maintainer: Ward9250)Development status:(Image: Build Status) (Image: Build status) (Image: codecov) (Image: )"
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
    "text": "Install NaturalSelection from the Julia REPL:julia> Pkg.add(\"NaturalSelection\")If you are interested in the cutting edge of the development, please check out the master branch to try new features before release."
},

{
    "location": "index.html#Contributing-and-Questions-1",
    "page": "Home",
    "title": "Contributing and Questions",
    "category": "section",
    "text": "We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features. Please go to the contributing section of the documentation for more information.If you have a question about contributing or using this package, you are encouraged to use the Bio category of the Julia discourse site."
},

{
    "location": "man/dNdS.html#",
    "page": "dNdS",
    "title": "dNdS",
    "category": "page",
    "text": "CurrentModule = NaturalSelection"
},

{
    "location": "man/dNdS.html#dN-/-dS-1",
    "page": "dNdS",
    "title": "dN / dS",
    "category": "section",
    "text": "NaturalSelection.jl provides several different methods of inferring the action of natural selection from coding sequences.Evolutionary pressures on proteins are often quantified by the ratio of substitution rates at non-synonymous and synonymous sites i.e. dN/dS.The dN/dS ratio was originally developed for application to distantly diverged sequences, the differences among which represent substitutions that have fixed along independent lineages.Nevertheless, the dN/dS measure is often applied to sequences sampled from a single population, the differences among which represent segregating polymorphisms. However, do be careful if this is what you are doing, as it has been demonstrated that dN/dS is not always suitable for such purposes (Sergey Kryazhimskiy & Joshua B. Plotkin, 2008)."
},

{
    "location": "man/dNdS.html#NaturalSelection.dNdS_NG86",
    "page": "dNdS",
    "title": "NaturalSelection.dNdS_NG86",
    "category": "Function",
    "text": "dNdS_NG86(x, y, addone::Bool = true, code::Int = 1)\n\nCompute dN and dS, using the Nei and Gojobori 1986 method.\n\nThe genetic code that is used, is defined according to the numbering of ncbi_trans_table. Code 1 is the standard genetic code.\n\nThis function requires two iterables x and y. If these iterables yield Codon{DNA} or Codon{RNA} type variables. Then it is assumed that x and y are iterables that yield a sequence of aligned codons. If the iterables produce DNA or RNA type variables, then it is assumed x and y iterables that conform to the behaviour of DNA or RNA sequences as defined in the BioSequences package. In this case, a new x and y that do have an element type of Codon{DNA} or Codon{RNA}.\n\nNG86 is a counting method of computing dN/dS and is typically safer to use on sequence data where codon usage, (esp. at 3rd position), is uniform, the sequences are not very divergent, and transition/transversion rates, are similar.\n\n\n\n"
},

{
    "location": "man/dNdS.html#NaturalSelection.S_N_NG86",
    "page": "dNdS",
    "title": "NaturalSelection.S_N_NG86",
    "category": "Function",
    "text": "S_N_NG86(codon::C, code::GeneticCode) where {C <: CDN}\n\nEnumerate the number of synonymous (S) and non-synonymous (N) sites in a codon, using the method used by Nei and Gojobori (1986).\n\nReturns a tuple where S is the first element and N is the second (S, N).\n\nEach site in a codon may be both partially synonymous and non-synonymous.\n\n\n\n"
},

{
    "location": "man/dNdS.html#NaturalSelection.DS_DN_NG86",
    "page": "dNdS",
    "title": "NaturalSelection.DS_DN_NG86",
    "category": "Function",
    "text": "DS_DN_NG86(x::C, y::C, code::GeneticCode) where C <: CDN\n\nCompute the number of synonymous (DS) and non-synonymous (DN) mutations between two codons, using the all paths method used by the Nei and Gojobori (1986).\n\n\n\n"
},

{
    "location": "man/dNdS.html#The-NG86-method-1",
    "page": "dNdS",
    "title": "The NG86 method",
    "category": "section",
    "text": "dNdS_NG86\nS_N_NG86\nDS_DN_NG86"
},

{
    "location": "man/tajimad.html#",
    "page": "Tajima's D",
    "title": "Tajima's D",
    "category": "page",
    "text": "CurrentModule = NaturalSelection"
},

{
    "location": "man/tajimad.html#NaturalSelection.tajimad-Tuple{AbstractFloat,Integer,Integer}",
    "page": "Tajima's D",
    "title": "NaturalSelection.tajimad",
    "category": "Method",
    "text": "tajimad(π::AbstractFloat, S::Integer, n::Integer)\n\nCompute Tajima's D from:\n\nπ: The average number of SNPs found in (n choose 2) pairwise comparisons of      a sample of sequences.\nS: The number of segregating sites in a sample of sequences.\nn: The number of sequences in your sample.\n\nExample\n\ntajimad(3.88888, 16, 10)\n\n\n\n"
},

{
    "location": "man/tajimad.html#NaturalSelection.tajimad-Tuple{Any}",
    "page": "Tajima's D",
    "title": "NaturalSelection.tajimad",
    "category": "Method",
    "text": "tajimad(seqs)\n\nCompute Tajima's D from a collection of BioSequences{DNAAlphabet{n}} (n = 2 or 4).\n\nThis will estimate the π, S, and n parameters from the sequences and use those parameters to estimate Tajima's D.\n\nExample\n\n\nsample = [dna\"ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA\",\n          dna\"AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA\",\n          dna\"AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA\",\n          dna\"AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA\",\n          dna\"AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA\",\n          dna\"AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA\",\n          dna\"AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA\",\n          dna\"AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA\",\n          dna\"AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA\",\n          dna\"AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA\"]\n\ntajimad(sample)\n\n\n\n"
},

{
    "location": "man/tajimad.html#Tajima's-D-1",
    "page": "Tajima's D",
    "title": "Tajima's D",
    "category": "section",
    "text": "Tajima's D is a population genetic test statistic created by and named after the Japanese researcher Fumio Tajima.Tajima's D is computed as the difference between two measures of genetic diversity: The mean number of pairwise differences and the number of segregating sites, each scaled so that they are expected to be the same in a neutrally evolving population of constant size.The purpose of the statistic is to distinguish between a DNA sequence evolving randomly (\"neutrally\") and one evolving under a non-random process. The non-random process might be directional or balancing selection, demographic expansion or contraction, genetic hitchhiking, or even introgression.tajimad(::AbstractFloat, ::Integer, ::Integer)\ntajimad(::Any)"
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
