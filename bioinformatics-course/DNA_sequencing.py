import copy


def string_composition(text, k):
    """This function returns all k-mers of a given text in lexicographic order."""
    composition_k = []  # to store composition
    for i in range(len(text) - k + 1):  # loop over the text
        composition_k.append(text[i:i+k])  # identify each k-mer
    return composition_k


def path_to_genome(path):
    """This function returns a text based on its genome path."""
    text = []  # to store text
    for i in range(len(path)-1):  # loop over the genome path, omit the last k-mer
        text.append(path[i][0])  # add first element of each k-mer to the text
    text.append(path[-1])  # join last k-mer to the text
    text = "".join(map(str, text))  # list to string
    return text


def cycle_to_genome(path):
    """This function returns a text based on its genome path."""
    text = []  # to store text
    for i in range(len(path) - 1):  # loop over the genome path, omit the last k-mer
        text.append(path[i][0])  # add first element of each k-mer to the text
    text = "".join(map(str, text))  # list to string
    return text


def prefix(pattern):
    """This is a function that returns the prefix of a pattern (last character is removed)."""
    return pattern[0:len(pattern) - 1]  # remove last character


def suffix(pattern):
    """This is a function that returns the suffix of a pattern (first character is removed)."""
    return pattern[1:len(pattern)]  # remove first character


def overlap_graph(pattern):
    """This function returns the overlap graph of a collection of k-mers (pattern)."""
    adjacency_list = dict()  # to store overlapping k-mers
    for i in range(len(pattern)):  # loop over the pattern
        suffices = []
        for j in range(len(pattern)):  # for each k-mer, loop over the pattern
            if suffix(pattern[i]) == prefix(pattern[j]):  # if the k-mer overlaps
                suffices.append(pattern[j])  # collect k-mer in list
        if len(suffices) > 0:  # if one or more overlapping k-mers are identified
            adjacency_list[pattern[i]] = suffices  # collect overlapping k-mers in dictionary
    return adjacency_list


def de_bruijn_graph(text, k):
    """This function returns the k-mer overlap graph of a given text taking into account gluing.
    If text is known!"""
    patterns = []  # to store k-mers
    for i in range(len(text) - k + 1):  # loop over text
        patterns.append(text[i:i + k])  # identify each k-mer
    adjacency_list = {prefix(kmer): [] for kmer in patterns}  # to store overlapping k-mers
    for i in range(len(patterns) - 1):  # loop over k-mers
        adjacency_list[prefix(patterns[i])].append(prefix(patterns[i + 1]))  # link prefix k-mer to next k-mer's prefix
    adjacency_list[prefix(patterns[-1])].append(suffix(patterns[-1]))  # link prefix last k-mer to last k-mer's suffix
    return adjacency_list


def de_bruijn(patterns):
    """This function returns a graph in which every k-mer in Patterns is isolated edge between its prefix and suffix.
    The graph results from ﻿gluing all nodes in dB with identical labels."""
    adjacency_list = dict()  # to store overlapping k-mers
    for kmer in patterns:  # add all possible unique nodes (k-1-mers) to dictionary
        adjacency_list[prefix(kmer)] = []
        adjacency_list[suffix(kmer)] = []
    for kmer in patterns:  # link prefix of each k-mer to it's suffix
        adjacency_list[prefix(kmer)].append(suffix(kmer))
    for node in list(adjacency_list):  # remove all unique nodes (k-1-mers) that have no overlapping node
        if not adjacency_list[node]:
            del adjacency_list[node]
    return adjacency_list


def cycle(node, graph):
    """This function returns a path by randomly walking in graph without visiting the same edge twice."""
    used_nodes = [node]  # to store nodes that have been used
    circle = []  # to store nodes that no longer have any connections (in cycle in reverse!)
    import copy
    unused_connections = copy.deepcopy(graph)  # to keep track of the connections that can still be used

    # to determine when the circle is finished:
    length_circle = 0
    for i in graph.keys():  # determine number of connections that should be made
        length_circle += len(graph[i])  # sum all connections

    # compute circle:
    while len(circle) < length_circle:  # as long as not all connections are used
        # as long as the current node still has unused connections
        while node in unused_connections and unused_connections[node]:
            old_node = node  # store current node
            node = unused_connections[node][0]  # update node to first possible connection of current node
            # keep track of the connections that can still be used
            unused_connections[old_node].remove(unused_connections[old_node][0])
            used_nodes.append(node)  # keep track of nodes that have been used
        circle.append(node)  # add node that has no connections left --> circle in reverse
        used_nodes.pop()  # remove this node from used_nodes
        node = used_nodes[-1]  # continue with previous node (if it still has connections left)
    circle.append(node)  # add last node to the circle
    circle.reverse()  # reverse the cycle
    return circle


def eulerian_cycle(graph):
    """This function returns a cycle by randomly walking in graph without visiting the same edge twice."""
    import random
    node = random.choice(list(graph))  # initialize node
    # node = "".join(map(str, [0]*(k-1)))
    circle = cycle(node, graph)  # compute path
    return circle


def eulerian_path(graph):
    """This function returns a path by randomly walking in graph without visiting the same edge twice."""
    graph_out = graph  # rename graph of outgoing connections
    graph_in = dict()  # initialize graph of incoming connections

    # create graph of incoming connections
    for node in graph:  # loop over each node in graph of outgoing connections
        for connection in range(len(graph[node])):  # loop over all outgoing connections in current node
            if graph[node][connection] in graph_in:  # if node is already in dictionary (for multiple connections)
                graph_in[graph[node][connection]].append(node)  # include all incoming connection
            if graph[node][connection] not in graph_in:  # if node is not yet in dictionary
                graph_in[graph[node][connection]] = []  # add node to dictionary
                graph_in[graph[node][connection]].append(node)  # include incoming connection

    # find starting node
    node_start = None
    for node in graph_out:  # loop over nodes in graph of outgoing connections
        if node not in graph_in:
            node_start = node
        elif len(graph_out[node]) > len(graph_in[node]):  # start at node with less outgoing than incoming connections
            node_start = node

    path = cycle(node_start, graph)  # compute path
    return path


def string_reconstruction(k, patterns):
    """This function assembles a genome from a given set of k-mers (patterns)"""
    graph = de_bruijn(patterns)  # construct the de Bruijn graph from the given pattern
    path = eulerian_path(graph)  # construct the Eulerian path from the computed graph
    text = path_to_genome(path)  # construct text from the computed path
    return text


def k_universal_circular_string(k):
    from itertools import product
    patterns = [''.join(p) for p in product('10', repeat=k)]
    patterns.sort()
    graph = de_bruijn(patterns)  # construct the de Bruijn graph from the given pattern
    cycle = eulerian_cycle(graph)  # construct the Eulerian path from the computed graph
    text = path_to_genome(cycle)  # construct text from the computed path
    text = text[0:len(text) - (k-1)]
    return text


def string_spelled_by_patterns(graph):
    """This function returns a text corresponding to a path in a de Bruijn graph."""
    path = eulerian_path(graph)  # construct the Eulerian path from the computed graph
    text = path_to_genome(path)  # construct text from the computed path
    return text


def string_spelled_by_gapped_patterns(patterns, k, d):
    """This function returns a sequence of (k,d)-mers corresponding to a path in a paired de Bruijn graph.
    Suffix((ai|bi)) = Prefix((ai+1|bi+1)) for 1 <= i <= n-1.
    len(Text) = k + d + k + n - 1. """
    patterns_split = list((line.strip().split('|') for line in patterns))  # split at |

    first_patterns = []
    second_patterns = []
    for i in range(len(patterns)):  # collect first and second pattern in separate list
        first_patterns.append(patterns_split[i][0])
        second_patterns.append(patterns_split[i][1])

    prefix_text = path_to_genome(first_patterns)  # compute genome from path of first patterns
    suffix_text = path_to_genome(second_patterns)  # compute genome from path of second patterns

    if prefix_text[k + d:] != suffix_text[:-k - d]:  # check if prefix and suffix of text overlap
        return "No such string"

    return prefix_text + suffix_text[-k - d:]


def de_bruijn_from_read_pairs(patterns):
    """This function returns a graph in which every k-mer in Patterns is isolated edge between its prefix and suffix.
    The graph results from ﻿gluing all nodes in dB with identical labels."""
    adjacency_list = dict()  # to store overlapping k-mers
    patterns_split = list((line.strip().split('|') for line in patterns))
    for i in range(len(patterns)):  # add all possible unique nodes (k-1-mer read-pairs) to dictionary
        adjacency_list[prefix(patterns_split[i][0]) + '|' + prefix(patterns_split[i][1])] = []
        adjacency_list[suffix(patterns_split[i][0]) + '|' + suffix(patterns_split[i][1])] = []
    for i in range(len(patterns)):  # link prefix of each (k-1-mer read pair) to it's suffix
        adjacency_list[prefix(patterns_split[i][0]) + '|' + prefix(patterns_split[i][1])].append(
            suffix(patterns_split[i][0]) + '|' + suffix(patterns_split[i][1]))
    for node in list(adjacency_list):  # remove all unique nodes (k-1-mers) that have no overlapping node
        if not adjacency_list[node]:
            del adjacency_list[node]
    return adjacency_list


def string_reconstruction_from_read_pairs(patterns, k, d):
    """This function reconstructs a text based on its paired composition of k,d-mers."""
    graph = de_bruijn_from_read_pairs(patterns)
    path = eulerian_path(graph)
    text = string_spelled_by_gapped_patterns(path, k, d)
    return text


def maximal_non_branching_paths(graph):
    graph_out = graph  # rename graph of outgoing connections
    graph_in = dict()  # initialize graph of incoming connections
    graph_tot = copy.deepcopy(graph_out)  # initialize graph of total connections

    # create graph of incoming connections & total connections
    for node in graph:  # loop over each node in graph of outgoing connections
        for connection in range(len(graph[node])):  # loop over all outgoing connections in current node
            if graph[node][connection] in graph_in:  # if node is already in dictionary (for multiple connections)
                graph_in[graph[node][connection]].append(node)  # include all incoming connection
            if graph[node][connection] not in graph_in:  # if node is not yet in dictionary
                graph_in[graph[node][connection]] = []  # add node to dictionary
                graph_in[graph[node][connection]].append(node)  # include incoming connection
            if graph[node][connection] in graph_tot:  # if node is already in dictionary (for multiple connections)
                graph_tot[graph[node][connection]].append(node)  # include all incoming connection
            if graph[node][connection] not in graph_tot:  # if node is not yet in dictionary
                graph_tot[graph[node][connection]] = []  # add node to dictionary
                graph_tot[graph[node][connection]].append(node)  # include incoming connection

    # find non-branching paths
    paths = []
    one_in_one_out = []
    used = []
    for node in graph_tot:
        if node in graph_out and node in graph_in:
            if len(graph_in[node]) == 1 and len(graph_in[node]) == len(graph_out[node]):
                one_in_one_out.append(node)
    for node in graph_tot:
        if node not in one_in_one_out:
            if node in graph_out and len(graph_out[node]) > 0:
                for node_w in graph_out[node]:
                    used.append(node)
                    used.append(node_w)
                    nbp = [node, node_w]
                    while node_w in one_in_one_out:
                        nbp.append(graph_out[node_w][0])
                        node_w = graph_out[node_w][0]
                        used.append(node_w)
                    paths.append(nbp)

    # find isolated cycles
    for node in graph_tot:
        if node not in used:
            used.append(node)
            isolated_cycle = [node]
            start = node
            flag = False
            while node in one_in_one_out:
                isolated_cycle.append(graph_out[node][0])
                node = graph_out[node][0]
                used.append(node)
                if node == start:
                    flag = True
                    break
            if flag:
                paths.append(isolated_cycle)

    # to reformat as graph
    paths_print = [0]*len(paths)
    for i in range(len(paths)):
        paths_print[i] = '->'.join(map(str, paths[i]))

    return paths


def contig_generation(patterns):
    """This function computes the contigs from a set of patterns. """
    contigs = []
    graph = de_bruijn(patterns)  # find the graph from the given set of patterns
    paths = maximal_non_branching_paths(graph)  # find the paths of all contigs
    for i in range(len(paths)):  # reformat to text
        contigs.append(path_to_genome(paths[i]))
    return contigs
