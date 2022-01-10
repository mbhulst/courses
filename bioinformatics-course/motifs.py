# NOTE: function names and values contain capitals because Stepik uses Pseudocode

from replication import HammingDistance
from replication import Neighbors
from replication import ApproximatePatternCount
from replication import NumberToPattern

import random


def GibbsSampler(Dna, k, t, N):
    """This function randomly selects the best motif in a given set of Texts (Dna)."""
    BestMotifs = RandomMotifs(Dna, k, t)
    for j in range(N):
        r = random.randint(0, t - 1)
        TempMotifs = BestMotifs
        Profile = ProfileWithPseudocounts(TempMotifs[:r] + TempMotifs[r + 1:])
        TempMotifs[r] = ProfileGeneratedString(Dna[r], Profile, k)
        if Score(TempMotifs) < Score(BestMotifs):
            BestMotifs = TempMotifs
    return BestMotifs


def ProfileGeneratedString(Text, Profile, k):
    """This function returns a randomly generated k-mer from Text whose probabilities are generated from Profile."""
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):  # loop over k-mers in the given Text
        # determine the probability that this k-mer is chosen based on the given profile
        probabilities[Text[i:i + k]] = Pr(Text[i:i + k], Profile)
    probabilities = Normalize(probabilities)  # normalize the probabilities
    return WeightedDie(probabilities)  # return the randomly chosen k-mer based on the computed probabilities


def WeightedDie(Probabilities):
    """This function returns a randomly chosen k-mer key with respect to the values in probabilities.
    Note: probabilities should sum to 1. """
    r = random.uniform(0, 1)  # choose a random decimal between 0 and 1
    probabilities_sum = 0
    for kmer in Probabilities:  # loop over k-mers in probability dictionary
        probabilities_sum += Probabilities[kmer]  # each iteration, increase by the probability of current k-mer
        if r < probabilities_sum:  # if the random value fall in this interval, return current k-mer
            return kmer


def Normalize(Probabilities):
    """This function rescales a collection of probabilities so that these probabilities sum to 1."""
    return {key: Probabilities[key] / sum(Probabilities.values()) for key in Probabilities}


def RandomizedMotifSearch(Dna, k, t):
    """This is a function that returns the most probably motif of k-mers in a list of t texts (Dna).
     It is based on a random input. """
    M = RandomMotifs(Dna, k, t)  # randomly choose a set of k-mer motifs from Dna
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)  # compute the profile corresponding to the chosen set of k-mer motifs
        M = Motifs(Profile, Dna)  # find the most probable collection of k-mers based on the profile
        if Score(M) < Score(BestMotifs):  # if the score is better than before, return motifs, otherwise repeat
            BestMotifs = M
        else:
            return BestMotifs


def RandomMotifs(Dna, k, t):
    """This function chooses a random k-mer from each of t texts (Dna)."""
    RandomMotifs = []
    for i in range(t):  # for each string of Text in Dna
        r = random.randint(0, len(Dna[0]) - k)
        RandomMotifs.append(Dna[i][r:r + k])  # randomly choose a k-mer
    return RandomMotifs


def Motifs(Profile, Dna):
    """This a function that returns the most probably collection of k-mers in each string of Dna.
    It is based by the Profile-most probable k-mers."""
    k = len(Profile["A"])
    Motifs = []
    for i in range(len(Dna)):  # for each string of Text in Dna
        Text = Dna[i]
        Motif = ProfileMostProbableKmer(Text, k, Profile)  # find the most probable k-mer based on the given profile
        Motifs.append(Motif)
    return Motifs


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    """This function will find the most probable motif of length k in a list of t texts (Dna).
    Taking into account pseudocounts. """
    BestMotifs = []
    for i in range(0, t):  # initialize the vector with the first k-mer from each string of Dna
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])  # length of each text in the list
    for i in range(n - k + 1):  # loop over the first string of Dna
        Motifs = []  # reset motifs
        Motifs.append(Dna[0][i:i + k])  # collect a motif with length k of first string of Dna
        for j in range(1, t):
            profile = ProfileWithPseudocounts(Motifs[0:j])  # determine profile of this motif
            # find the most probable motif in the other strings of Dna that matches with the motif in the first string
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, profile))  # collect these in motifs
        if Score(Motifs) < Score(BestMotifs):  # if the score is better, this will be best motif
            BestMotifs = Motifs
    return BestMotifs


def ProfileWithPseudocounts(Motifs):
    """This function provides a matrix of counts of each nucleotide in a list of k-mers divided by
    the number of k-mers provided (t), the profile of the motif. Taking into account pseudocounts. """
    t = len(Motifs)  # number of k-mers in the list
    k = len(Motifs[0])  # length of the k-mers
    profile = {symbol: [0] * k for symbol in "ACGT"}  # initialize the dictionary
    count = CountWithPseudocounts(Motifs)  # determine the matrix of counts
    for symbol in "ACGT":  # for each nucleotide
        for i in range(k):  # for each k-mer index
            profile[symbol][i] = count[symbol][i] / (t + 4)  # divide the character by the total number of k-mers (t)
    return profile


def CountWithPseudocounts(Motifs):
    """ This function provides a matrix of counts of each nucleotide in a list of k-mers, the count of the motif.
     Taking into account Pseudocounts. """
    t = len(Motifs)  # number of k-mers in the list
    k = len(Motifs[0])  # length of the k-mers
    count = {symbol: [1] * k for symbol in "ACGT"}  # initialize the dictionary
    for i in range(t):  # for each row
        for j in range(k):  # for each column
            symbol = Motifs[i][j]  # for each k-mer, loop over each nucleotide
            count[symbol][j] += 1  # count the nucleotides in the right dictionary key
    return count


def MedianString(Dna, k):
    """This function will find the median pattern of lenght k in a given set of Dna sequences."""
    distance = 1E10
    for i in range(4 ** k):  # loop over all possible Patterns with length k
        Pattern = NumberToPattern(i, k)
        if distance > DistanceBetweenPatternAndStrings(Pattern, Dna):
            distance = DistanceBetweenPatternAndStrings(Pattern, Dna)
            Median = Pattern
    return Median


def DistanceBetweenPatternAndStrings(Pattern, Dna):
    """This function calculates the sum of distances between Pattern and each Text in Dna."""
    t = len(Dna)  # number to texts in Dna
    n = len(Dna[0])  # length of texts in Dna
    k = len(Pattern)  # length of the pattern
    distance = 0
    for i in range(t):  # loop over each text in Dna
        Hammingdistance = 1E10
        for j in range(n - k + 1):  # loop over each possible pattern in each text in Dna
            window = Dna[i][j: j+k]  # define window
            if Hammingdistance > HammingDistance(Pattern, window):  # if current pattern is closer to given pattern
                Hammingdistance = HammingDistance(Pattern, window)  # this is the best pattern
        distance = distance + Hammingdistance  # summarize best hamming distances
    return distance


def MotifEnumeration(Dna, k, d):
    """This function returns the (k,d)-motif that appears in every string from Dna with at most d mismatches."""
    Patterns = []
    n = len(Dna[0])  # length of each text in the list
    for i in range(n - k + 1):  # loop over the first text in Dna
        Text = Dna[0][i: i+k]  # define window
        if d == 0:
            Pattern = Text
            count = 0
            for l, dna in enumerate(Dna):  # loop over each text in the list
                if ApproximatePatternCount(dna, Pattern, d) >= 1:  # if the pattern is at least once in the current list
                    count = count + 1  # store in count
                if count == len(Dna):  # if lists contain a pattern, store in patterns
                    Patterns.append(Pattern)
        else:
            for j, Pattern in enumerate(Neighbors(Text, d)):  # loop over each pattern in the neighborhood of the window
                count = 0
                for l, dna in enumerate(Dna):  # loop over each text in the list
                    if ApproximatePatternCount(dna, Pattern, d) >= 1:  # if the pattern is at least once in the current list
                        count = count + 1  # store in count
                    if count == len(Dna):  # if lists contain a pattern, store in patterns
                        Patterns.append(Pattern)
    Patterns = list(dict.fromkeys(Patterns))
    return Patterns


def GreedyMotifSearch(Dna, k, t):
    """This function will find the most probable motif of length k in a list of t texts (Dna)"""
    BestMotifs = []
    for i in range(0, t):  # initialize the vector with the first k-mer from each string of Dna
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])  # length of each text in the list
    for i in range(n - k + 1):  # loop over the first string of Dna
        Motifs = []  # reset motifs
        Motifs.append(Dna[0][i:i + k])  # collect a motif with length k of first string of Dna
        for j in range(1, t):
            profile = Profile(Motifs[0:j])  # determine profile of this motif
            # find the most probable motif in the other strings of Dna that matches with the motif in the first string
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, profile))  # collect these in motifs
        if Score(Motifs) < Score(BestMotifs):  # if the score is better, this will be best motif
            BestMotifs = Motifs
    return BestMotifs


def ProfileMostProbableKmer(Text, k, Profile):
    """This function finds the most probable k-mer in a text with a given profile matrix."""
    score = -1
    for i in range(len(Text) - k + 1):  # Loop over the text
        Pattern = Text[i: i+k]  # Define the window
        p = Pr(Pattern, Profile)
        if p > score:
            kmer = Pattern
            score = p
    return kmer


def Pr(Text, Profile):
    """This function calculates the probability that a profile matrix will produce a given string."""
    p = 1
    for i in range(len(Text)):  # for each column (index of profile)
        p *= Profile[Text[i]][i]
    return p


def Score(Motifs):
    """This function provides the score of a given set of motifs, based on the consensus sequence."""
    t = len(Motifs)  # number of k-mers in the list
    k = len(Motifs[0])  # length of the k-mers
    score = 0
    consensus = Consensus(Motifs)
    for i in range(t):  # for each row (each k-mer)
        for j in range(k):  # for each column (each index of motif)
            if Motifs[i][j] != consensus[j]:  # if the character is not similar to the consensus -> +1
                score += 1
    return score


def Consensus(Motifs):
    """This function provides the consensus sequence of a list of k-mers. """
    k = len(Motifs[0])  # length of the k-mers
    count = Count(Motifs)  # determine the matrix of counts
    consensus = ""  # initialize the consensus sequence
    for j in range(k):  # for each k-mer index
        m = 0  # rest the count
        frequent_symbol = ""  # reset the symbol
        for symbol in "ACGT":  # for each nucleotide
            if count[symbol][j] > m:  # find nucleotide with highest count in current k-mer index
                m = count[symbol][j]  #
                frequent_symbol = symbol
        consensus += frequent_symbol  # add that nucleotide to the consensus sequence
    return consensus


def Profile(Motifs):
    """This function provides a matrix of counts of each nucleotide in a list of k-mers divided by
    the number of k-mers provided (t), the profile of the motif. """
    t = len(Motifs)  # number of k-mers in the list
    k = len(Motifs[0])  # length of the k-mers
    profile = {symbol: [0] * k for symbol in "ACGT"}  # initialize the dictionary
    count = Count(Motifs)  # determine the matrix of counts
    for symbol in "ACGT":  # for each nucleotide
        for i in range(k):  # for each k-mer index
            profile[symbol][i] = count[symbol][i]/t  # divide the character by the total number of k-mers (t)
    return profile


def Count(Motifs):
    """ This function provides a matrix of counts of each nucleotide in a list of k-mers, the count of the motif. """
    t = len(Motifs)  # number of k-mers in the list
    k = len(Motifs[0])  # length of the k-mers
    count = {symbol: [0] * k for symbol in "ACGT"}  # initialize the dictionary
    for i in range(t):  # for each row
        for j in range(k):  # for each column
            symbol = Motifs[i][j]  # for each k-mer, loop over each nucleotide
            count[symbol][j] += 1  # count the nucleotides in the right dictionary key
    return count
