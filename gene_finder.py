# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T':
        return 'A'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    dna = dna[::-1]
    comp = ''
    for a in range(len(dna)):
        comp = comp + get_complement(dna[a])
    return comp

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    i = 0
    while i<len(dna):
        if (dna[i:i+3] == 'TGA') or (dna[i:i+3] == 'TAG') or (dna[i:i+3]== 'TAA'):
            break
        i = i+3
    return dna[0:i]

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    i = 0
    count = 0
    orf_oneframe = []
    while i<len(dna):
        if dna[i:i+3] == 'ATG':
            orf_oneframe.append(rest_of_ORF(dna[i:]))
            i = len(orf_oneframe[count])+i
            count = count +1
        i = i+3
    return orf_oneframe

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    orf_threeframes = []
    orf_threeframes.extend(find_all_ORFs_oneframe(dna[0:]))
    orf_threeframes.extend(find_all_ORFs_oneframe(dna[1:]))
    orf_threeframes.extend(find_all_ORFs_oneframe(dna[2:]))
    return orf_threeframes

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    orf_both = []
    orf_both.extend(find_all_ORFs(dna))
    orf_both.extend(find_all_ORFs(get_reverse_complement(dna)))
    return orf_both

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orf_list = find_all_ORFs_both_strands(dna)
    longest = ''
    for i in orf_list:
        if len(i)>len(longest):
            longest = i
    return longest

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest = 0
    for i in range(num_trials):
        dna = shuffle_string(dna)
        current = len(longest_ORF(dna))
        if current > longest:
            longest = current
    return longest

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    i = 0
    aminos = ''
    while i < len(dna):
        if len(dna)-i >= 3:
            dna_fragment = dna[i:i+3]
            aminos = aminos + aa_table[dna_fragment]
        i += 3
    return aminos

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    dna_list = find_all_ORFs_both_strands(dna)
    amino_list = []
    threhold = longest_ORF_noncoding(dna,1500)
    for i in range(len(dna_list)):
        if len(dna_list[i])>= threhold:
            amino = coding_strand_to_AA(dna_list[i])
            amino_list.append(amino)
    return amino_list
    #return coding_strand_to_AA(longest_ORF(dna))

if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(coding_strand_to_AA,globals(),verbose = True)
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
    #doctest.testmod()
