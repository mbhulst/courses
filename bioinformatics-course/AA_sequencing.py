from replication import ReverseComplement
import copy


def DNA_to_RNA(text):
    text_RNA = text.replace('T', 'U')
    return text_RNA


def RNA_to_DNA(text):
    text_DNA = text.replace('U', 'T')
    return text_DNA


def protein_translation(pattern, genetic_code):
    """This function translates a given RNA string (pattern) into a amino acid string."""
    peptide = ""
    i = 0
    while i < len(pattern) - 2:
        codon = pattern[i:i + 3]
        peptide += genetic_code[codon]
        i += 3
    return peptide


def encoding_problem(text, peptide, genetic_code):
    """This function returns all possible sequences in a given DNA string (text) that could encode a given peptide. """
    answer = []  # initialize answer
    DNA_seq = text  # rename DNA sequence
    DNA_seq_RC = ReverseComplement(text)  # take reverse complement of DNA sequence

    # define the 6 reading frames of the RNA sequence
    DNA_seqs = []
    DNA_seqs.append(DNA_seq)
    DNA_seqs.append(DNA_seq[1:])
    DNA_seqs.append(DNA_seq[2:])
    DNA_seqs.append(DNA_seq_RC)
    DNA_seqs.append(DNA_seq_RC[1:])
    DNA_seqs.append(DNA_seq_RC[2:])

    for frame in range(len(DNA_seqs)):  # for each reading frame
        RNA_seq = DNA_to_RNA(DNA_seqs[frame])  # compute RNA sequence
        AA_seq = protein_translation(RNA_seq, genetic_code)  # compute amino acid sequence
        i = 0
        while i < len(AA_seq):  # loop over de amino acid sequence of this frame to find all sequences of peptide
            i = AA_seq.find(peptide, i)  # find the given amino acid sequence (peptide)
            if i == -1:  # stop if it is not found
                break
            if frame in range(0, 3):  # when it is found in the forward frames, store sequence
                answer.append(DNA_seqs[frame][i * 3: i * 3 + len(peptide) * 3])
            elif frame in range(3, 6):  # when it is found in the reverse complement (RC) frames, store RC sequence
                answer.append(ReverseComplement(DNA_seqs[frame][i * 3: i * 3 + len(peptide) * 3]))
            i += len(peptide)
    return answer


def generating_linear_spectrum(peptide, integer_mass_table):
    """This function generates a linear spectrum of subpeptides for a given peptide sequence. """
    prefix_mass = [0]*(len(peptide) + 1)  # to store the masses of the prefixes of each amino acid in the sequence
    for i in range(1, len(peptide) + 1):  # loop over the peptide
        prefix_mass[i] = prefix_mass[i-1] + integer_mass_table[peptide[i-1]]  # calculate mass of prefix of current AA
    linear_spectrum = [0]  # to store the linear spectrum
    for i in range(len(peptide) + 1):  # loop over peptide
        for j in range(i+2, len(peptide) + 2):  # loop over all amino acids after current AA
            linear_spectrum.append(prefix_mass[j-1] - prefix_mass[i])  # calculate mass of the peptide [i:j-1]
    return sorted(linear_spectrum)


def generating_cyclic_spectrum(peptide, integer_mass_table):
    """This function generates a cyclic spectrum of subpeptides for a given peptide sequence. """
    prefix_mass = [0] * (len(peptide) + 1)  # to store the masses of the prefixes of each amino acid in the sequence
    for i in range(1, len(peptide) + 1):  # loop over the peptide
        prefix_mass[i] = prefix_mass[i - 1] + integer_mass_table[peptide[i - 1]]  # calculate mass of prefix of AA
    peptide_mass = prefix_mass[-1]
    cyclic_spectrum = [0]  # to store the linear spectrum
    for i in range(len(peptide) + 1):  # loop over peptide
        for j in range(i + 2, len(peptide) + 2):  # loop over all amino acids after current AA
            cyclic_spectrum.append(prefix_mass[j - 1] - prefix_mass[i])  # calculate mass of the peptide [i:j-1]
            if i > 0 and j - 1 < len(peptide):  # also add subpeptides that wrap around the end of peptide
                cyclic_spectrum.append(peptide_mass - (prefix_mass[j - 1] - prefix_mass[i]))
    return sorted(cyclic_spectrum)


def counting_peptides_with_given_mass(mass, integer_mass_table):
    """This function returns the number of linear peptides having a given mass."""
    del integer_mass_table["L"]  # remove L because I/L are considered the same
    del integer_mass_table["Q"]  # remove Q because K/Q are considered the same

    routes = [0]*(mass + 1)
    routes[0] = 1
    for i in range(mass + 1):
        for j in integer_mass_table:
            if i + integer_mass_table[j] < mass + 1:
                routes[i + integer_mass_table[j]] += routes[i]
    count = routes[-1]
    return count


def counting_subpeptides(mass):
    """This function returns the number of subpeptides of a linear peptide."""
    count = 1 + mass
    for i in range(mass):
        count += i
    return count


def mass_peptide(peptide, integer_mass_table):
    """This function returns the mass of a given peptide."""
    mass = 0
    for i in range(len(peptide)):  # loop over each amino acid in peptide and add its mass
        mass += integer_mass_table[peptide[i]]
    return mass


def expand_peptides(peptides):
    """This function expands a list of peptides with each possible amino acid. """
    peptides_expanded = []
    for i in range(len(peptides)):  # loop over each peptide in the list
        for j in integer_mass_table:  # expand with each possible amino acid
            peptides_expanded.append(peptides[i] + j)
    return peptides_expanded


def check_spectrum(small_list, large_list):
    """This function checks if list B contains all items of list A (including duplicates)."""
    B_check = copy.deepcopy(large_list)  # the list should not be changed in the rest of the code!
    for i in range(len(small_list)):  # loop over all items in small list
        if small_list[i] in B_check:  # if the item is found in large list, remove it from large list
            B_check.remove(small_list[i])
        else:  # if the item is not found --> does not contain all items of small list
            return False
    return True  # if all items can be removed --> does contain all items of small list


def cyclopeptide_sequencing(spectrum):
    """This function returns a cyclic peptide whose theoretical spectrum matches the given experimental spectrum."""
    del integer_mass_table["L"]  # remove L because I/L are considered the same
    del integer_mass_table["Q"]  # remove Q because K/Q are considered the same
    candidate_peptides = [""]  # to store the peptides that have to be checked
    final_peptides = []  # to collect the peptides that match the given spectrum

    while candidate_peptides:  # as long as some candidate peptides remain
        candidate_peptides = expand_peptides(candidate_peptides, integer_mass_table)  # expand from current list

        # NOTE: loop in reverse, otherwise list is changed during loop due to deletion!
        for peptide in reversed(candidate_peptides):  # loop over each candidate peptide
            mass = mass_peptide(peptide, integer_mass_table)  # determine mass of current peptide
            spectrum_check_linear = generating_linear_spectrum(peptide, integer_mass_table)  # compute linear spectrum
            spectrum_check_cyclic = generating_cyclic_spectrum(peptide, integer_mass_table)  # compute cyclic spectrum

            # if the last mass of the spectrum is found the peptide could be a match
            if mass == spectrum[-1]:  # is the mass the same?
                if spectrum_check_cyclic == spectrum and peptide not in final_peptides:  # is the spectrum the same?
                    final_peptides.append(peptide)  # store the peptide
                candidate_peptides.remove(peptide)  # if so, remove the peptide from candidates

            # check if all fragments occur in the spectrum
            # (NOTE: use linear spectrum, because it is a spectrum of fragments!)
            elif check_spectrum(spectrum_check_linear, spectrum) is False:  # are all fragments in the spectrum?
                candidate_peptides.remove(peptide)  # if not, remove the peptide from candidates

    # formatting for stepik
    final_peptides_formatted = []
    for i in range(len(final_peptides)):
        masses = []
        for j in range(len(final_peptides[i])):
            masses.append(mass_peptide(final_peptides[i][j], integer_mass_table))
        masses = '-'.join(map(str, masses))  # if the result needs to be split by an arrow
        final_peptides_formatted.append(masses)
    return final_peptides_formatted


def score_cyclic(peptide, spectrum):
    """This function returns the score of an experimental spectrum against the theoretical spectrum
    of a cyclic peptide."""
    score = 0  # initialize score
    spectrum_exp = copy.deepcopy(spectrum)  # make copy because we are removing from the list
    spectrum_theo = generating_cyclic_spectrum(peptide, integer_mass_table)  # generate theoretical spectrum
    for i in range(len(spectrum_theo)):  # loop over each fragment in theoretical spectrum
        if spectrum_theo[i] in spectrum_exp:  # if it is found in the experimental spectrum
            spectrum_exp.remove(spectrum_theo[i])  # remove the fragment
            score += 1  # and count the score
    return score


def score_linear(peptide, spectrum):
    """This function returns the score of an experimental spectrum against the theoretical spectrum
    of a linear peptide."""
    score = 0  # initialize score
    spectrum_exp = copy.deepcopy(spectrum)  # make copy because we are removing from the list
    spectrum_theo = generating_linear_spectrum(peptide, integer_mass_table)  # generate theoretical spectrum
    for i in range(len(spectrum_theo)):  # loop over each fragment in theoretical spectrum
        if spectrum_theo[i] in spectrum_exp:  # if it is found in the experimental spectrum
            spectrum_exp.remove(spectrum_theo[i])  # remove the fragment
            score += 1  # and count the score
    return score


def trim(leaderboard, spectrum, N):
    """This function returns the top N highest-scoring peptides in leaderboard with respect to the spectrum."""
    linear_scores = [0]*len(leaderboard)

    for i in range(len(leaderboard)):  # loop over the leaderboard
        peptide = leaderboard[i]  # pick a peptide
        linear_scores[i] = score_linear(peptide, spectrum)  # determine its (LINEAR) score

    # sort leaderboard according to scores & reverse
    leaderboard = [x for (y, x) in sorted(zip(linear_scores, leaderboard))]
    leaderboard = leaderboard[::-1]

    # sort & reverse linear scores
    linear_scores = sorted(linear_scores)
    linear_scores = linear_scores[::-1]

    for i in range(N, len(leaderboard)):  # loop over leaderboard
        if linear_scores[i] < linear_scores[N]:
            leaderboard = leaderboard[:N]  # remove everything with a lower score than peptide N
            return leaderboard
    return leaderboard


def leaderboard_cyclopeptide_sequencing(spectrum, N):
    """This function returns a cyclic peptide having maximum score against an experimental spectrum."""
    del integer_mass_table["L"]  # remove L because I/L are considered the same
    del integer_mass_table["Q"]  # remove Q because K/Q are considered the same

    leaderboard = [""]  # to store the N highest scoring linear peptides including ties
    leader_peptide = ""  # to collect the peptide that matches the given spectrum
    storage = []

    while leaderboard:  # while leaderboard is non-empty
        leaderboard = expand_peptides(leaderboard)  # expand from current list

        # NOTE: loop in reverse, otherwise list is changed during loop due to deletion!
        for peptide in reversed(leaderboard):  # loop over each candidate peptide
            mass = mass_peptide(peptide, integer_mass_table)  # determine mass of current peptide

            # if score_cyclic(peptide, spectrum) == 83:
            #     storage.append(peptide)

            # if the mass is same as highest mass in spectrum
            if mass == spectrum[-1]:  # is the mass the same?
                if score_linear(peptide, spectrum) > score_linear(leader_peptide, spectrum):  # is the score higher?
                    leader_peptide = peptide  # replace the peptide
                    print(peptide, score_linear(peptide, spectrum), score_cyclic(peptide, spectrum))
            # if the mass is not the same as highest mass in spectrum
            elif mass > spectrum[-1]:
                leaderboard.remove(peptide)  # remove the peptide from candidates
        leaderboard = trim(leaderboard, spectrum, N)

    masses = []
    for i in range(len(leader_peptide)):
        masses.append(mass_peptide(leader_peptide[i], integer_mass_table))
    leader_peptide_formatted = '-'.join(map(str, masses))  # if the result needs to be split by an arrow
    print(leader_peptide)
    print(score_cyclic(leader_peptide, spectrum))
    # print(storage)
    return leader_peptide_formatted


# construct dictionary of amino acid integer masses
with open('integer_mass_table.txt', 'r') as file:  # if data is given as graph
    integer_mass_table = {line.strip().split(' ')[0]: line.strip().split(' ')[1] for line in file}
integer_mass_table = {k: int(v) for k, v in integer_mass_table.items()}

# construct dictionary of genetic code
with open('RNA_codon_table.txt', 'r') as file:  # if data is given as graph
    genetic_code = dict(line.strip().split(' ') for line in file)
# invert_genetic_code = {dest: [k for k, v in genetic_code.items() if v == dest] for dest in set(genetic_code.values())}

with open('data.txt', 'r') as data:  # if the data is split by spaces
    x = data.read().split()
x = [int(y) for y in x]

# spectrum = x
# N = 1000

# answer = leaderboard_cyclopeptide_sequencing(spectrum, N)
# print(answer)
