
###########################################################################

def mkt(list_of_list_of_list_of_seq):
    ''' (list) -> table, float
    McDonald-Kreitman test.
    >>> mkt([[['AAAGGGCCC', 'AGAGGGCGC'],['AACGGGCCG','ACAGGGCCC']], [['ATGGGA', 'AGAGGA'],['ATGCGA']]])
                   Polymorphism     Divergence
    Synonymous     3                0
    Non-synonymous 4                2
    Alpha = 1.0
    '''
    a = []; b = []; c = []; d = []; e =[]; g = []; ps = []; pn = []; ds = []; dn = []
    for i in range(len(list_of_list_of_list_of_seq)):
        a.append(poly_ratio(list_of_list_of_list_of_seq[i][0]))
        b.append(poly_ratio(list_of_list_of_list_of_seq[i][1]))
        c.append(div_ratio(list_of_list_of_list_of_seq[i][0], list_of_list_of_list_of_seq[i][1]))
    for i in range(len(a)):
        ps.append(a[i][0])
        ps.append(b[i][0])
        d.append(a[i][0]+b[i][0]+1)
    for i in range(len(a)):
        pn.append(a[i][1])
        pn.append(b[i][1])
        e.append(a[i][1]+b[i][1])
    for i in range(len(c)):
        ds.append(c[i][0])
        dn.append(c[i][1])
    if len(list_of_list_of_list_of_seq) == 1:
        if sum(dn) and sum(ps) != 0:
            alpha = 1 - ((sum(ds)*sum(pn))/(sum(dn)*sum(ps)))
        else:
            alpha = 'NULL'
    elif len(list_of_list_of_list_of_seq) > 1:
        for i in range(len(d)):
            g.append(e[i]/d[i])
        if sum(dn) != 0:
            alpha = 1 - (sum(ds)/sum(dn))*(sum(g)/len(list_of_list_of_list_of_seq))
        else:
            alpha = 'NULL'
    print('               Polymorphism     Divergence' + '\n' +
          'Synonymous     ' + str(sum(ps)) + str((17-len(str(sum(ps))))*' ') + str(sum(ds)) + '\n' +
          'Non-synonymous ' + str(sum(pn)) + str((17-len(str(sum(pn))))*' ') + str(sum(dn)) + '\n' +
          'Alpha = ' + str(alpha))


###########################################################################

def poly_ratio(list_of_seq):
    ''' (list) -> list
    Returns the total number of synonymous and non-synonymous mutations between
    a list of aligned coding sequences. The first value of the output list is
    the total number of synonymous mutations and the second value of the output
    list is the total number of non-synonymous mutations.
    >>> poly_ratio(['ATGCACAAAGGG', 'ATGCTCAGAGGG', 'ATGCACAAGGGG'])
    [1, 2]
    >>> poly_ratio(['AAAGGGCCC', 'AGAGGGCGC', 'AACGGGCCG','ACAGGGCCC'])
    [1, 4]
    '''
    a = codon_diff(list_of_seq); p = [0, 0]
    for i in range(len(a)):
        if len(a[i]) > 1:
            p[0] += multi_short_path(a[i])[0]
            p[1] += multi_short_path(a[i])[1]
    return p

###########################################################################

def div_ratio(list_of_seq_1, list_of_seq_2):
    ''' (list, list) -> list
    Returns the total number of synonymous and non-synonymous divergent changes
    between two lists of aligned sequences. list_of_seq_1 belongs to species_1 and
    list_of_seq_2 belongs to species_2. The first value of the output list is the
    total number of synonymous divergent changes and the second value of the output
    list is the total number of non-synonymous divergent changes.
    >>> div_ratio(['ATGCACAAAGGG', 'ATGCTCAGAGGG', 'ATGCACAAGGGG'], ['ATGCACGGAGGG', 'ATGTTAAGGGGG'])
    [1, 0]
    '''
    a = codon_diff(list_of_seq_1); b = codon_diff(list_of_seq_2); d = [0, 0]
    for i in range(len(a)):
        if len(list(set(a[i]) & set(b[i]))) == 0:
            x = []
            for j in range(len(b[i])):
                for k in range(len(a[i])):
                    x.append(short_path(b[i][j], a[i][k]))
            d[0] += helper(x)[0][0]
            d[1] += helper(x)[0][1]
    return d

###########################################################################

def codon_diff(list_of_seq):
    ''' (list) -> dict
    Given a list of aligned sequences, it returns a dictionary with keys as codon
    positions and values as the variants of that codon.
    >>> codon_diff(['AAAGGGCCC', 'AGAGGGCGC', 'AACGGGCCG','ACAGGGCCC'])
    {0: ['AAA', 'AGA', 'AAC', 'ACA'], 1: ['GGG'], 2: ['CCC', 'CGC', 'CCG']}
    '''
    d = {}
    for i in range(len(list_of_seq[0])):
        if i % 3 == 0:
            d[i//3] = [list_of_seq[0][i:(i+3)]]
    for i in range(len(list_of_seq)):
        for j in range(len(list_of_seq[0])):
            if list_of_seq[0][j] != list_of_seq[i][j]:
                if j % 3 == 0:
                    a = list_of_seq[i][j:(j+3)]
                    if a not in d[j//3]:
                        d[j//3].append(a)
                elif j % 3 == 1:
                    a = list_of_seq[i][(j-1):(j+2)]
                    if a not in d[j//3]:
                        d[j//3].append(a)
                elif j % 3 == 2:
                    a = list_of_seq[i][(j-2):(j+1)]
                    if a not in d[j//3]:
                        d[j//3].append(a)
    return d

###########################################################################

def multi_short_path(list_of_codon):
    ''' (list) -> list
    Returns the shortest mutation pathway between a list of codons as a list of two
    integers. The first integer of the list is the amount of synonymous mutations and the
    second is the number of non-synonymous mutations required for the most probable
    mutation pathway between codons. The pathway with the least overall amount of
    nucleotide substitutions (synonymous + non-synonymous) and the least amount of
    non-synonymous mutations is considered as the most probable. Codons with incomplete
    sequences are excluded from the final computation. The computation is made using
    Kruskal's alogrithm.
    >>> multi_short_path(['AAA', 'AGA', 'AAG'])
    [1, 1]
    >>> multi_short_path(['AAA', 'AGA', 'AAG', 'CCC'])
    [3, 2]
    '''
    mat = []; p = [0,0]; g = []; h = []
    for i in range(len(list_of_codon)):
        for j in range(len(list_of_codon)):
            mat.append(short_path(list_of_codon[i], list_of_codon[j]))
            g.append(list_of_codon[i])
            g.append(list_of_codon[j])
    for k in range(len(mat)):
        a = (get_rank(mat[k]), g[k*2], g[(k*2)+1])
        h.append(a)
    graph = {'vertices':list_of_codon, 'edges':set(h)}
    t = kruskal(graph)
    f = list(t)
    for m in range(len(f)):
        n = get_lst(f[m][0])
        p[0] += n[0]
        p[1] += n[1]
    return p

###########################################################################

def short_path(codon_1, codon_2):
    ''' (str, str) -> list
    Returns the shortest mutation pathway between two codons as a list of two
    integers. The first number of the list is the amount of synonymous mutations
    and the second is the number of non-synonymous mutations required for the most
    probable mutation pathway between two codons. The pathway with the least
    overall amount of nucleotide substitutions (synonymous + non-synonymous) and
    the least amount of non-synonymous mutations is considered as the most probable.
    >>> short_path('GTG', 'AAA')
    [1, 2]
    >>> short_path('ATG', 'CTA')
    [1, 1]
    '''
    lis = []; las =[]
    if codon_1 == codon_2:
        return [0, 0]
    for i in range(len(codon_1)):
        if codon_1[i] != codon_2[i]:
            a = ''
            a = codon_1[:i] + codon_2[i] + codon_1[(i+1):]
            lis.append(a)
    for i in range(len(lis)):
        for j in range(len(codon_2)):
            if codon_2[j] != lis[i][j]:
                a = ''
                a = lis[i][:j] + codon_2[j] + lis[i][(j+1):]
                las.append(a)
    if lis == []:
        return [0,0]
    elif codon_2 in lis:
        if nuc_to_aa(codon_1) == nuc_to_aa(codon_2):
            return [1, 0]
        else:
            return [0, 1]
    elif codon_2 in las:
        if (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) == nuc_to_aa(codon_2)):
            return [2,0]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) != nuc_to_aa(codon_2)):
            return [1, 1]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[0]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[1]) == nuc_to_aa(codon_2)):
            return [1, 1]
        else:
            return [0, 2]
    else:
        if (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) == nuc_to_aa(las[0]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) == nuc_to_aa(las[1]) == nuc_to_aa(codon_2)):
            return [3, 0]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) == nuc_to_aa(las[2]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) == nuc_to_aa(las[3]) == nuc_to_aa(codon_2)):
            return [3, 0]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[2]) == nuc_to_aa(las[4]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[2]) == nuc_to_aa(las[5]) == nuc_to_aa(codon_2)):
            return [3, 0]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) == nuc_to_aa(las[0]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) == nuc_to_aa(las[1]) != nuc_to_aa(codon_2)):
            return [2, 1]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) == nuc_to_aa(las[2]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) == nuc_to_aa(las[3]) != nuc_to_aa(codon_2)):
            return [2, 1]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[2]) == nuc_to_aa(las[4]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[2]) == nuc_to_aa(las[5]) != nuc_to_aa(codon_2)):
            return [2, 1]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) != nuc_to_aa(las[0]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) != nuc_to_aa(las[1]) == nuc_to_aa(codon_2)):
            return [2, 1]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) != nuc_to_aa(las[2]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) != nuc_to_aa(las[3]) == nuc_to_aa(codon_2)):
            return [2, 1]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[2]) != nuc_to_aa(las[4]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[2]) != nuc_to_aa(las[5]) == nuc_to_aa(codon_2)):
            return [2, 1]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[0]) == nuc_to_aa(las[0]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[0]) == nuc_to_aa(las[1]) == nuc_to_aa(codon_2)):
            return [2, 1]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[1]) == nuc_to_aa(las[2]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[1]) == nuc_to_aa(las[3]) == nuc_to_aa(codon_2)):
            return [2, 1]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[2]) == nuc_to_aa(las[4]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[2]) == nuc_to_aa(las[5]) == nuc_to_aa(codon_2)):
            return [2, 1]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[0]) != nuc_to_aa(las[0]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[0]) != nuc_to_aa(las[1]) == nuc_to_aa(codon_2)):
            return [1, 2]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[1]) != nuc_to_aa(las[2]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[1]) != nuc_to_aa(las[3]) == nuc_to_aa(codon_2)):
            return [1, 2]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[2]) != nuc_to_aa(las[4]) == nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[2]) != nuc_to_aa(las[5]) == nuc_to_aa(codon_2)):
            return [1, 2]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[0]) == nuc_to_aa(las[0]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[0]) == nuc_to_aa(las[1]) != nuc_to_aa(codon_2)):
            return [1, 2]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[1]) == nuc_to_aa(las[2]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[1]) == nuc_to_aa(las[3]) != nuc_to_aa(codon_2)):
            return [1, 2]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[2]) == nuc_to_aa(las[4]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[2]) == nuc_to_aa(las[5]) != nuc_to_aa(codon_2)):
            return [1, 2]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) != nuc_to_aa(las[0]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[0]) != nuc_to_aa(las[1]) != nuc_to_aa(codon_2)):
            return [1, 2]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) != nuc_to_aa(las[2]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[1]) != nuc_to_aa(las[3]) != nuc_to_aa(codon_2)):
            return [1, 2]
        elif (nuc_to_aa(codon_1) == nuc_to_aa(lis[2]) != nuc_to_aa(las[4]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) == nuc_to_aa(lis[2]) != nuc_to_aa(las[5]) != nuc_to_aa(codon_2)):
            return [1, 2]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[0]) != nuc_to_aa(las[0]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[0]) != nuc_to_aa(las[1]) != nuc_to_aa(codon_2)):
            return [0, 3]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[1]) != nuc_to_aa(las[2]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[1]) != nuc_to_aa(las[3]) != nuc_to_aa(codon_2)):
            return [0, 3]
        elif (nuc_to_aa(codon_1) != nuc_to_aa(lis[2]) != nuc_to_aa(las[4]) != nuc_to_aa(codon_2)) or (nuc_to_aa(codon_1) != nuc_to_aa(lis[2]) != nuc_to_aa(las[5]) != nuc_to_aa(codon_2)):
            return [0, 3]

###########################################################################

def helper(list_of_list):
    ''' (list) -> (list)
    Helper function for the div_ratio() function. Sorts a list of lists of two
    integers. The inner lists represent the number synonymous and non-synonymous
    mutations, where the first integer of an inner list is the amount of synonymous
    mutations and the second is the number of non-synonymous mutations required
    for the most probable mutation pathway between two codons. The sorting is done
    by the following hieararchy: [1, 0], [0, 1], [2, 0], [1, 1], [0, 2], [3, 0],
    [2, 1], [1, 2], [0, 3]. This hierarchy represents all possible substitutions
    between codons regarding the amount of synonymous and non-synonynmous mutations
    required. The first list of the hierarchy is the [1, 0] list which represents
    one synonymous and zero non-synonymous mutations. The hierarchy is continued
    giving priority to the list with the least amount of overall nucleotide
    substitutions (synonynmous + non-synonymous) and the least amount of
    non-synonymous substitutions.
    >>> helper([[1, 1], [0, 1]])
    [[0, 1], [1, 1]]
    >>> helper([[0, 1], [1, 0], [2, 0], [0, 2], [1, 0], [1, 1]])
    [[1, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]]
    '''
    a = []
    for i in range(len(list_of_list)):
        if list_of_list[i] == [1, 0]:
            a.append(list_of_list[i])
    for i in range(len(list_of_list)):
        if list_of_list[i] == [0, 1]:
            a.append(list_of_list[i])
    for i in range(len(list_of_list)):
        if list_of_list[i] == [2, 0]:
            a.append(list_of_list[i])
    for i in range(len(list_of_list)):
        if list_of_list[i] == [1, 1]:
            a.append(list_of_list[i])
    for i in range(len(list_of_list)):
        if list_of_list[i] == [0, 2]:
            a.append(list_of_list[i])
    for i in range(len(list_of_list)):
        if list_of_list[i] == [3, 0]:
            a.append(list_of_list[i])
    for i in range(len(list_of_list)):
        if list_of_list[i] == [2, 1]:
            a.append(list_of_list[i])
    for i in range(len(list_of_list)):
        if list_of_list[i] == [1, 2]:
            a.append(list_of_list[i])
    for i in range(len(list_of_list)):
        if list_of_list[i] == [0, 3]:
            a.append(list_of_list[i])
    return a

###########################################################################

def get_rank(lst):
    ''' (list) -> int
    Helper function for the multi_short_path() function. Ranks all possible
    nucleotide substitution combinations between two codons. The first number of
    the lst represent the number of synonymous substitutions while the second
    represents the number of non-synonymous substitutions. The ranking is done by
    the following hieararchy: [1, 0], [0, 1], [2, 0], [1, 1], [0, 2], [3, 0],
    [2, 1], [1, 2], [0, 3].
    >>> get_rank([0, 1])
    2
    >>> get_rank([0, 2])
    5
    '''
    rank = 0
    if lst == [1, 0]:
        rank = 1
    elif lst == [0, 1]:
        rank = 2
    elif lst == [2, 0]:
        rank = 3
    elif lst == [1, 1]:
        rank = 4
    elif lst == [0, 2]:
        rank = 5
    elif lst == [3, 0]:
        rank = 6
    elif lst == [2, 1]:
        rank = 7
    elif lst == [1, 2]:
        rank = 8
    elif lst == [0, 3]:
        rank = 9
    return rank

###########################################################################

def get_lst(rank):
    ''' (int) -> list
    Helper function for the multi_short_path() function. A reverse function to the
    get_rank()function.
    >>> get_lst(2)
    [0, 1]
    >>> get_lst(5)
    [0, 2]
    '''
    lst = [0, 0]
    if rank == 1:
        lst = [1, 0]
    elif rank == 2:
        lst = [0, 1]
    elif rank == 3:
        lst = [2, 0]
    elif rank == 4:
        lst = [1, 1]
    elif rank == 5:
        lst = [0, 2]
    elif rank == 6:
        lst = [3, 0]
    elif rank == 7:
        lst = [2, 1]
    elif rank == 8:
        lst = [1, 2]
    elif rank == 9:
        lst = [0, 3]
    return lst

###########################################################################

def nuc_to_aa(seq):
    ''' (str) -> str
    Returns a string of amino acids coded by the sequence.
    >>> nuc_to_aa('ATGGCCATG')
    'MAM'
    '''
    aa = ''
    for i in range(len(seq)):
        if i%3 == 0:
            aa += str((get_aa(seq[i:(i+3)])))
    return aa

###########################################################################

def get_aa(codon):
    ''' (str) -> str
    Returns the aminoacid coded by the codon.
    >>> get_aa('AUG')
    'M'
    >>> get_aa('CAA')
    'Q'
    '''
    if codon in ['UUU', 'UUC', 'TTT', 'TTC']:
        return 'F'
    elif codon in ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']:
        return 'L'
    elif codon in ['AUU', 'AUC', 'AUA', 'ATT', 'ATC', 'ATA']:
        return 'I'
    elif codon in ['AUG', 'ATG']:
        return 'M'
    elif codon in ['GUU', 'GUC', 'GUA', 'GUG', 'GTT', 'GTC', 'GTA', 'GTG']:
        return 'V'
    elif codon in ['UCU', 'UCC', 'UCA', 'UCG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGU', 'AGC', 'AGT']:
        return 'S'
    elif codon in ['CCU', 'CCC', 'CCA', 'CCG', 'CCT']:
        return 'P'
    elif codon in ['ACU', 'ACC', 'ACA', 'ACG', 'ACT']:
        return 'T'
    elif codon in ['GCU', 'GCC', 'GCA', 'GCG', 'GCT']:
        return 'A'
    elif codon in ['UAU', 'UAC', 'TAT', 'TAC']:
        return 'Y'
    elif codon in ['CAU', 'CAC', 'CAT']:
        return 'H'
    elif codon in ['CAA', 'CAG']:
        return 'Q'
    elif codon in ['AAU', 'AAC', 'AAT']:
        return 'N'
    elif codon in ['AAA', 'AAG']:
        return 'K'
    elif codon in ['GAU', 'GAC', 'GAT']:
        return 'D'
    elif codon in ['GAA', 'GAG']:
        return 'E'
    elif codon in ['UGU', 'UGC', 'TGT', 'TGC']:
        return 'C'
    elif codon in ['UGG', 'TGG']:
        return 'W'
    elif codon in ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG', 'CGT']:
        return 'R'
    elif codon in ['GGU', 'GGC', 'GGA', 'GGG', 'GGT']:
        return 'G'
    elif codon in ['UAA', 'UAG', 'UGA', 'TAA', 'TAG', 'TGA']:
        return '*'

###########################################################################


def pathways_a(codon_1, codon_2):
    ''' (str, str) -> print
    Returns all possible direct mutation pathways between two codons in
    aminoacid form.
    >>> pathways_a('AAA', 'AGA')
    K -> R
    >>> pathways_a('ATG', 'CTA')
    M -> L -> L
    M -> I -> L
    '''
    lis = []; las = []
    for i in range(len(codon_1)):
        if codon_1[i] != codon_2[i]:
            a = ''
            a = codon_1[:i] + codon_2[i] + codon_1[(i+1):]
            lis.append(a)
    for i in range(len(lis)):
        for j in range(len(codon_2)):
            if codon_2[j] != lis[i][j]:
                a = ''
                a = lis[i][:j] + codon_2[j] + lis[i][(j+1):]
                las.append(a)
    if lis == []:
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(codon_2))
    elif codon_2 in lis:
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(codon_2))
    elif codon_2 in las:
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(lis[0]) + ' -> ' + nuc_to_aa(codon_2))
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(lis[1]) + ' -> ' + nuc_to_aa(codon_2))
    else:
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(lis[0]) + ' -> ' + nuc_to_aa(las[0]) + ' -> ' + nuc_to_aa(codon_2))
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(lis[0]) + ' -> ' + nuc_to_aa(las[1]) + ' -> ' + nuc_to_aa(codon_2))
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(lis[1]) + ' -> ' + nuc_to_aa(las[2]) + ' -> ' + nuc_to_aa(codon_2))
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(lis[1]) + ' -> ' + nuc_to_aa(las[3]) + ' -> ' + nuc_to_aa(codon_2))
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(lis[2]) + ' -> ' + nuc_to_aa(las[4]) + ' -> ' + nuc_to_aa(codon_2))
        print(nuc_to_aa(codon_1) + ' -> ' + nuc_to_aa(lis[2]) + ' -> ' + nuc_to_aa(las[5]) + ' -> ' + nuc_to_aa(codon_2))

###########################################################################

def pathways_n(codon_1, codon_2):
    ''' (str, str) -> print
    Returns all possible direct mutation pathways between two codons in
    nucleotide form.
    >>> pathways_n('AAA', 'AGA')
    AAA -> AGA
    >>> pathways_n('ATG', 'CTA')
    ATG -> CTG -> CTA
    ATG -> ATA -> CTA
    '''
    lis = []; las = []
    for i in range(len(codon_1)):
        if codon_1[i] != codon_2[i]:
            a = ''
            a = codon_1[:i] + codon_2[i] + codon_1[(i+1):]
            lis.append(a)
    for i in range(len(lis)):
        for j in range(len(codon_2)):
            if codon_2[j] != lis[i][j]:
                a = ''
                a = lis[i][:j] + codon_2[j] + lis[i][(j+1):]
                las.append(a)
    if lis == []:
        print(codon_1 + ' -> ' + codon_2)
    elif codon_2 in lis:
        print(codon_1 + ' -> ' + codon_2)
    elif codon_2 in las:
        print(codon_1 + ' -> ' + lis[0] + ' -> ' + codon_2)
        print(codon_1 + ' -> ' + lis[1] + ' -> ' + codon_2)
    else:
        print(codon_1 + ' -> ' + lis[0] + ' -> ' + las[0] + ' -> ' + codon_2)
        print(codon_1 + ' -> ' + lis[0] + ' -> ' + las[1] + ' -> ' + codon_2)
        print(codon_1 + ' -> ' + lis[1] + ' -> ' + las[2] + ' -> ' + codon_2)
        print(codon_1 + ' -> ' + lis[1] + ' -> ' + las[3] + ' -> ' + codon_2)
        print(codon_1 + ' -> ' + lis[2] + ' -> ' + las[4] + ' -> ' + codon_2)
        print(codon_1 + ' -> ' + lis[2] + ' -> ' + las[5] + ' -> ' + codon_2)

###########################################################################

# Kruskal's algorithm - adopted from 'https://github.com/israelst/Algorithms-Book--Python/blob/master/5-Greedy-algorithms/kruskal.py'


parent = {}; rank = {}

def make_set(vertice):
    parent[vertice] = vertice
    rank[vertice] = 0

def find(vertice):
    if parent[vertice] != vertice:
        parent[vertice] = find(parent[vertice])
    return parent[vertice]

def union(vertice1, vertice2):
    root1 = find(vertice1)
    root2 = find(vertice2)
    if root1 != root2:
        if rank[root1] > rank[root2]:
            parent[root2] = root1
        else:
            parent[root1] = root2
            if rank[root1] == rank[root2]: rank[root2] += 1

def kruskal(graph):
    for vertice in graph['vertices']:
        make_set(vertice)
    minimum_spanning_tree = set()
    edges = list(graph['edges'])
    edges.sort()
    for edge in edges:
        weight, vertice1, vertice2 = edge
        if find(vertice1) != find(vertice2):
            union(vertice1, vertice2)
            minimum_spanning_tree.add(edge)
    return minimum_spanning_tree

###########################################################################

if __name__ == '__main__':
    import doctest
    doctest.testmod()
