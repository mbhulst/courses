# NOTE: function names and values contain capitals because Stepik uses Pseudocode

def ApproximatePatternCount(Text, Pattern, d):
    """This function that returns the # of occurrences of patterns with distance <= d to given pattern in given text."""
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        window = Text[i:i + len(Pattern)]  # define window
        distance = HammingDistance(Pattern, window)  # find distance between window and given pattern
        if distance <= d:
            count += 1  # count # time pattern with distance <= d occurs
    return count


def ApproximatePatternMatching(Text, Pattern, d):
    """This function returns the starting positions in text of each pattern with distance <= d to the given pattern."""
    positions = []
    for i in range(len(Text) - len(Pattern) + 1):
        window = Text[i:i + len(Pattern)]  # define window
        distance = HammingDistance(Pattern, window)  # find distance between window and given pattern
        if distance <= d:
            positions.append(i)  # save position if distance <= d
    return positions


def HammingDistance(p, q):
    """This function returns the distance between to patterns p and q."""
    distance = 0
    n = len(p)
    for i in range(n):
        if p[i] != q[i]:  # where p is different than q
            distance += 1
    return distance


def MinimumSkew(Genome):
    """This is a function that returns all integers i minimizing Skew[i] for genome"""
    positions = []
    n = len(Genome)
    skew = SkewArray(Genome)
    min_value = min(skew)
    for i in range(n):
        if skew[i] == min_value:
            positions.append(i)
    return positions


def SkewArray(Genome):
    """This is a function that returns the skew array of a genome."""
    skew = [0] * (len(Genome) + 1)
    n = len(Genome)
    skew[0] = 0

    for i in range(n):
        if Genome[i] == "A":
            skew[i + 1] = skew[i]
        if Genome[i] == "C":
            skew[i + 1] = skew[i] - 1
        if Genome[i] == "G":
            skew[i + 1] = skew[i] + 1
        if Genome[i] == "T":
            skew[i + 1] = skew[i]
    return skew


def FasterSymbolArray(Genome, symbol):
    """This is a function that returns the symbol array of a text (genome) corresponding to a symbol."""

    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n // 2]

    # look at the first half of genome to compute first array value
    array[0] = PatternCount(ExtendedGenome[0:n // 2], symbol)

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i - 1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i - 1] == symbol:
            array[i] = array[i] - 1
        if ExtendedGenome[i + (n // 2) - 1] == symbol:
            array[i] = array[i] + 1
    return array


def SymbolArray(Genome, symbol):
    """This is a function that returns the symbol array of a text (genome) corresponding to a symbol."""
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n // 2]  # Paste half of the genome at the end
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i + (n // 2)], symbol)  # Count # of symbols in window
    return array


def BetterClumpFinding(Genome, k, L, t):
    """This is a function that finds all distinct k-mers (patterns) forming (L,t)-clumps in a given genome.
     Here t is the amount of times a k-mer should occur in a part of the genome with length L."""
    FrequentPatterns = []
    Clump = [0] * (4 ** k)
    Text = Genome[0:L]  # Define first window
    FrequencyArray = ComputingFrequencies(Text, k)  # Compute frequency array of first window
    for i in range(0, 4 ** k):
        if FrequencyArray[i] >= t:
            Clump[i] = 1  # Remember k-mers with frequency < t in first window
    for i in range(1, len(Genome) - L + 1):  # Loop over remaining windows
        FirstPattern = Genome[i - 1:i - 1 + k]  # Define first k-mer (pattern) of previous window
        index = PatternToNumber(FirstPattern)  # Compute index of that pattern
        FrequencyArray[index] = FrequencyArray[index] - 1  # Remove that pattern from frequency array
        LastPattern = Genome[i + L - k:i + L]  # Define last k-mer (pattern) of current window
        index = PatternToNumber(LastPattern)  # Compute index of that pattern
        FrequencyArray[index] = FrequencyArray[index] + 1  # Add that pattern to frequency array
        if FrequencyArray[index] >= t:
            Clump[index] = 1  # Remember k-mers with frequency < t in current window
    for i in range(0, 4 ** k):
        if Clump[i] == 1:  # Find index of all patterns with occurrence >= t
            Pattern = NumberToPattern(i, k)  # Find patterns belonging to those indexes
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
    # Return patterns (k-mers) that occur more than t time within L characters


def ClumpFinding(Genome, k, L, t):
    """This is a function that finds all distinct k-mers (patterns) forming (L,t)-clumps in a given genome.
     Here t is the amount of times a k-mer should occur in a part of the genome with length L.
     NOTE: this function is very slow!"""
    FrequentPatterns = []
    Clump = [0] * (4 ** k)
    for i in range(0, len(Genome) - L + 1):  # Loop over genome with length L
        Text = Genome[i:i + L]  # Define window
        FrequencyArray = ComputingFrequencies(Text, k)  # Compute frequency array of window
        for index in range(0, len(FrequencyArray) - k + 1):  # Loop over frequency array
            if FrequencyArray[index] >= t:
                Clump[index] = 1  # Update index numbers to 1 when occurrence >= t
    for i in range(0, 4 ** k):
        if Clump[i] == 1:  # Find index of all patterns with occurrence >= t
            Pattern = NumberToPattern(i, k)  # Find patterns belonging to those indexes
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
    # Return patterns (k-mers) that occur more than t time within L characters


def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
    """This is a function to determine the most frequent words (patterns) of length 'k' in a text
    accounting for mismatches of a given distance d and reverse complements."""
    FrequentPatterns = []
    FrequencyArray_MM = ComputingFrequenciesWithMismatches(Text, k, d)  # Compute the frequency array
    FrequencyArray_RC = ComputingFrequenciesWithMismatches(ReverseComplement(Text), k, d)
    FrequencyArray = [sum(x) for x in zip(FrequencyArray_MM, FrequencyArray_RC)]
    maxCount = max(FrequencyArray)  # Find the maximum # of occurrences of a pattern in the Text
    for i in range(0, len(FrequencyArray) - k + 1):  # Loop over the frequency array
        if FrequencyArray[i] == maxCount:
            # Find the index (number) of pattern with the maximum # occurrences
            Pattern = NumberToPattern(i, k)  # Find the pattern belonging to that index (number)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns  # Return the patterns that occur most frequently


def FrequentWordsWithMismatches(Text, k, d):
    """This is a function to determine the most frequent words (patterns) of length 'k' in a text
    accounting for mismatches of a given distance d."""
    FrequentPatterns = []
    FrequencyArray = ComputingFrequenciesWithMismatches(Text, k, d)  # Compute the frequency array
    maxCount = max(FrequencyArray)  # Find the maximum # of occurrences of a pattern in the Text
    for i in range(0, len(FrequencyArray) - k + 1):  # Loop over the frequency array
        if FrequencyArray[i] == maxCount:
            # Find the index (number) of pattern with the maximum # occurrences
            Pattern = NumberToPattern(i, k)  # Find the pattern belonging to that index (number)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns  # Return the patterns that occur most frequently


def ComputingFrequenciesWithMismatches(Text, k, d):
    """This is a function to generate a frequency array of a given text for a given k-mer
    accounting for mismatches of a given distance d."""
    FrequencyArray = [0] * (4 ** k)
    for i in range(0, len(Text) - k + 1):  # Loop over the text
        Pattern = Text[i:i + k]  # Define the window
        Neighborhood = Neighbors(Pattern, d)  # Find neighborhood of the pattern
        if type(Neighborhood) == str:
            j = PatternToNumber(Neighborhood)
            FrequencyArray[j] += 1
        else:
            for x in range(0, len(Neighborhood)):
                ApproximatePattern = Neighborhood[x]  # for each neighbor
                j = PatternToNumber(ApproximatePattern)
                # Find the index (number) that belongs to the pattern of the current window
                FrequencyArray[j] += 1  # count the occurrence of the pattern
    return FrequencyArray


def Neighbors(Pattern, d):
    """This is a recursive function that finds the neighbors of a given DNA sequence (pattern)
    with a given distance d. """
    if d == 0:  # if distance = 0, return Pattern
        return Pattern
    if len(Pattern) == 1:  # if Pattern has one character, return all possible characters
        return ["A", "C", "G", "T"]
    Neighborhood = []
    SuffixNeighbors = Neighbors(Suffix(Pattern), d)  # find the neighbors of the suffix of pattern with distance d
    for i in range(0, len(SuffixNeighbors)):
        Text = SuffixNeighbors[i]
        if HammingDistance(Suffix(Pattern), Text) < d:  # if current pattern
            Neighborhood.append("A" + Text)
            Neighborhood.append("C" + Text)
            Neighborhood.append("G" + Text)
            Neighborhood.append("T" + Text)
        else:
            Neighborhood.append(Pattern[0] + Text)
    return Neighborhood


def Suffix(Pattern):
    """This is a function that returns the suffix of a pattern (first character is removed)."""
    return Pattern[1:len(Pattern)]


def FasterFrequentWords(Text, k):
    """This is a function to determine the most frequent words (patterns) of length 'k' in a text.
    It is faster than the function FrequentWords for small 'k'."""
    FrequentPatterns = []
    FrequencyArray = ComputingFrequencies(Text, k)  # Compute the frequency array
    maxCount = max(FrequencyArray)  # Find the maximum # of occurrences of a pattern in the Text
    for i in range(0, len(FrequencyArray) - k + 1):  # Loop over the frequency array
        if FrequencyArray[i] == maxCount:
            # Find the index (number) of pattern with the maximum # occurrences
            Pattern = NumberToPattern(i, k)  # Find the pattern belonging to that index (number)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns  # Return the patterns that occur most frequently


def ComputingFrequencies(Text, k):
    """This is a function to generate a frequency array of a given text for a given k-mer.
    This is a faster method than using the function PatternCount."""
    FrequencyArray = [0] * (4 ** k)
    for i in range(0, len(Text) - k + 1):  # Loop over the text
        Pattern = Text[i:i + k]  # Define the window
        j = PatternToNumber(Pattern)
        # Find the index (number) that belongs to the pattern of the current window
        FrequencyArray[j] += 1  # Count the occurrence of the pattern
    return FrequencyArray


def PatternCount(Text, Pattern):
    """This is a function that counts the number of times a pattern occurs in a given text. """
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count = count + 1
    return count


def FrequentWords(Text, k):
    """This is a function that finds the patterns (words) of length k (k-mers) that occur most in a given text."""
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
        # add each key to words whose corresponding frequency value is equal to m
    return words


def FrequencyMap(Text, k):
    """This is a function that finds all unique patterns (words) of length k (k-mers).
     It also calculates how often each unique pattern occurs in a given text."""
    freq = {}
    # dictionary of all unique k-mers and the frequency of their occurrence in the text
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i + k]
        if pattern in freq:
            freq[pattern] += 1
        else:
            freq[pattern] = 1
    return freq


def NumberToPattern(index, k):
    """This is a recursive function to find the pattern (DNA sequence) for a given index number (of a DNA sequence)."""
    if k == 1:
        return NumberToSymbol(index)  # When k = 1, find last nucleotide and stop
    else:
        prefixIndex = int(index / 4)  # Quotient
        r = index % 4  # Remainder
        symbol = NumberToSymbol(r)  # Find the symbol of the remainder
        PrefixPattern = NumberToPattern(prefixIndex, k - 1)
        return PrefixPattern + symbol  # Return the symbols for each step and keep finding the next symbols


def PatternToNumber(Pattern):
    """This is a recursive function to find the index number for a given pattern (DNA sequence)."""
    if Pattern == "":
        return 0
    else:
        symbol = Pattern[-1]  # What is the last symbol?
        prefix = Pattern[:-1]  # Sequence without the last symbol
    return 4 * PatternToNumber(prefix) + SymbolToNumber(symbol)


def PatternMatching(Pattern, Genome):
    """This is a function to find where a given pattern (DNA sequence) occurs in a given text (genome).
    It returns the index (position) at which the patterns occur."""
    positions = []
    for i in range(len(Genome) - len(Pattern) + 1):  # Loop over the text (genome)
        window = Genome[i:i + len(Pattern)]
        if window == Pattern:  # Check if the current window is equal to the given pattern
            positions.append(i)  # If so, return the index (position) of the current window
    return positions


def SymbolToNumber(symbol):
    """This is a function to find the symbol (DNA nucleotide) for a given index number (of a DNA sequence)."""
    if symbol == "A":
        return 0
    elif symbol == "C":
        return 1
    elif symbol == "G":
        return 2
    elif symbol == "T":
        return 3


def NumberToSymbol(number):
    """This is a function to find the index number for a given DNA nucleotide."""
    if number == 0:
        return "A"
    elif number == 1:
        return "C"
    elif number == 2:
        return "G"
    elif number == 3:
        return "T"


def ReverseComplement(Pattern):
    """This is a function to find how the reverse complement of a given DNA sequence (pattern)."""
    Rev_Pattern = Reverse(Pattern)  # First compile the reverse of the sequence
    Rev_Comp_Pattern = Complement(Rev_Pattern)  # Then compile the complement of the reverse sequence
    return Rev_Comp_Pattern


def Reverse(Pattern):
    """This is a function to find how the reverse of a given DNA sequence (pattern)."""
    rev = ""
    for char in range(len(Pattern)):  # Loop over all characters (nucleotides) in the pattern (DNA sequence)
        rev = Pattern[char] + rev  # Reverse the pattern (DNA sequence)
    return rev


def Complement(Pattern):
    """This is a function to find how the complement of a given DNA sequence (pattern)."""
    comp = ""
    for char in range(len(Pattern)):  # Loop over all characters (nucleotides) in the pattern (DNA sequence)
        if Pattern[char] == 'A':  # Replace each character (nucleotide) with its complement
            comp = comp + 'T'
        elif Pattern[char] == 'C':
            comp = comp + 'G'
        elif Pattern[char] == 'G':
            comp = comp + 'C'
        elif Pattern[char] == 'T':
            comp = comp + 'A'
    return comp
