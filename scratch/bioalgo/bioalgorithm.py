import argparse
import sys
import math
import random
import bisect
from pprint import pprint


VERSION = "0.0.1"

__author__ = 'vador'


def all_strings(alphabet, length):
    """Find the list of all strings of 'alphabet' of length 'length'"""

    c = [[]]
    for i in range(length):
        c = [[x] + y for x in alphabet for y in c]
    return c


def dstmatch(pattern, refpattern):
    """ Give distance of pattern from refpattern """
    if len(pattern) != len(refpattern):
        return -1
    nbdiff = 0
    for i in range(len(pattern)):
        if pattern[i] != refpattern[i]:
            nbdiff += 1
    return nbdiff


def appear_in_dna(dnastr, kmer, maxdist):
    res = False
    l = len(kmer)

    for i in range(len(dnastr) - l + 1):
        if dstmatch(kmer, dnastr[i:i + l]) <= maxdist:
            return True
    return res


def get_composition_iter(dnastr, k):
    l = len(dnastr)

    for i in range(l - k + 1):
        yield dnastr[i:i + k]


def get_composition(dnastr, k):
    i_compo = get_composition_iter(dnastr, k)

    return sorted(i_compo)


def getarguments():
    parser = argparse.ArgumentParser(description='Compute most frequent kmer in string')
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)

    args = parser.parse_args()
    content = args.infile.readlines()
    return content


def parse_extra_dataset(string_list, endwithcolumn=True):
    inp = "Input"
    outp = "Output"
    if endwithcolumn:
        inp = 'Input:'
        outp = 'Output:'
    if string_list[0].rstrip() != inp:
        return -1, -1
    string_list = map(lambda x: x.rstrip(), string_list)

    outpos = string_list.index(outp)
    input = string_list[1:outpos]
    output = string_list[outpos + 1:]
    return input, output


def appear(listDna, length, dist):
    """ return list of kmers appearing in each strand"""
    allkmerslist = all_strings(['A', 'C', 'G', 'T'], length)

    #print allkmerslist
    res = []
    for k in allkmerslist:
        #k = ['A', 'A', 'A', 'T', 'C']
        cntAppear = 0
        for i in listDna:
            cntAppear += appear_in_dna(i, k, dist)

        if cntAppear == len(listDna):
            res.append("".join(k))

    print " ".join(res)
    print len(res)


def profile_matrix(DnaArray):
    """ return an array containing profile matrix : order is ATCG"""
    m = len(DnaArray)
    n = len(DnaArray[0])
    result = [[0 for i in range(n)] for j in range(4)]
    for i in range(m):
        for j in range(n):
            tmpNuc = DnaArray[i][j]
            if tmpNuc == 'A':
                result[0][j] += 1
            elif tmpNuc == 'T':
                result[1][j] += 1
            elif tmpNuc == 'C':
                result[2][j] += 1
            elif tmpNuc == 'G':
                result[3][j] += 1

    for i in range(4):
        for j in range(n):
            result[i][j] = result[i][j] * 1.0 / m

    return result


def profile_matrix_pseudocount(DnaArray):
    """ return an array containing profile matrix : order is ATCG"""
    m = len(DnaArray)
    n = len(DnaArray[0])
    result = [[1 for i in range(n)] for j in range(4)]
    for i in range(m):
        for j in range(n):
            tmp_nuc = DnaArray[i][j]
            if tmp_nuc == 'A':
                result[0][j] += 1
            elif tmp_nuc == 'T':
                result[1][j] += 1
            elif tmp_nuc == 'C':
                result[2][j] += 1
            elif tmp_nuc == 'G':
                result[3][j] += 1

    for i in range(4):
        for j in range(n):
            result[i][j] = result[i][j] * 1.0 / (m + 1)

    return result


def profile_matrix_dict(profile_matrix_arr):
    """ Convert an array containing profile matrix into a dict """
    result = dict()

    result['A'] = profile_matrix_arr[0]
    result['T'] = profile_matrix_arr[1]
    result['C'] = profile_matrix_arr[2]
    result['G'] = profile_matrix_arr[3]
    return result


def entropy(profile_matrix):
    """ return entropy from a profile matrix """
    m = len(profile_matrix)
    n = len(profile_matrix[0])
    res = 0
    for j in range(n):
        colRes = 0
        for i in range(m):
            tmp = profile_matrix[i][j]
            if tmp > 0:
                colRes += -tmp * math.log(tmp, 2)
        res += colRes
    return res


def score(motif):
    """ return an array containing profile matrix : order is ATCG"""
    m = len(motif)
    n = len(motif[0])
    result = [[0 for i in range(4)] for j in range(n)]
    # Build count for each nucleotide at each position
    for i in range(m):
        for j in range(n):
            tmpNuc = motif[i][j]
            if tmpNuc == 'A':
                result[j][0] += 1
            elif tmpNuc == 'T':
                result[j][1] += 1
            elif tmpNuc == 'C':
                result[j][2] += 1
            elif tmpNuc == 'G':
                result[j][3] += 1
    count = 0
    # discard most frequent, and sum
    for i in result:
        tmpPopular = sorted(i, reverse=True)
        count += sum(tmpPopular[1:])
    return count


def dist_kmer_strand(kmer, strand):
    """ return kmer in strand most approaching input kmer """
    l = len(kmer)
    res = len(kmer)
    strandlength = len(strand)
    for i in range(strandlength - l):
        tmpDna = strand[i:i + l]
        tmpDst = dstmatch(kmer, tmpDna)
        if tmpDst < res:
            res = tmpDst
    return res


def dist_kmer_dnalist(kmer, listDna):
    """ return distance of a kmer from a set of dnaStrands
    (sum of distances of most approaching Kmer in each strand) """
    res = 0
    for i in listDna:
        res += dist_kmer_strand(kmer, i)
    return res


def find_bestcandidate_bruteforce(listDna, k):
    """ returns kmer minimizing distance for all strands in listDna"""
    candidates = all_strings(['A', 'C', 'G', 'T'], k)
    #listDna = ["AAATTGACGCAT","GACGACCACGTT","CGTCAGCGCCTG","GCTGAGCACCGG","AGTACGGGACAG"]
    minscore = len(listDna) * k
    bestcandidate = []
    for i in candidates:
        tmpScore = dist_kmer_dnalist(i, listDna)
        if tmpScore < minscore:
            minscore = tmpScore
            bestcandidate = i
        elif tmpScore == minscore:
            bestcandidate.append(i)

    print map("".join, bestcandidate)
    print minscore
    return bestcandidate


def probability_from_profile(kmer, profileMatrixDict):
    """ return probability for a given kmer within the profile distribution """
    res = 1
    for i in range(len(kmer)):
        res *= profileMatrixDict[kmer[i]][i]
    return res


def fill_profile_dict(profilelist):
    """ fill a dict with profile values as strings
    ( each column is a nucleotide) """
    profileDict = dict()

    nuclOrder = profilelist[0].split(" ")

    tmpArr = []
    for i in range(1, len(profilelist)):
        tmpVal = profilelist[i].split(" ")
        tmpArr.append(tmpVal)

    for i in range(len(nuclOrder)):
        profileDict[nuclOrder[i]] = map(lambda x: float(x[i]), tmpArr)
    return profileDict


def get_best_probability_for_profile(dnaStrand, length, profileDict):
    """ returns the length-mer in dnaStrand with maximum probability according to profile"""
    res = 0
    candidate = dnaStrand[0:length]
    for i in range(len(dnaStrand) - length + 1):
        tmpStr = dnaStrand[i:i + length]
        tmpPr = probability_from_profile(tmpStr, profileDict)
        if tmpPr > res:
            res = tmpPr
            candidate = tmpStr

    return candidate, res


def greedy_motif_search(dnaStr, k, t):
    candidateMotif = [dnaStr[i][0:k] for i in range(len(dnaStr))]
    candidateScore = k * t
    for m0 in range(len(dnaStr[0]) - k + 1):
        tmpMotif = [dnaStr[0][m0:m0 + k]]
        for i in range(1, t):
            tmpProfil = profile_matrix_dict(profile_matrix(tmpMotif))
            candidatei, candidateprobai = get_best_probability_for_profile(dnaStr[i], k, tmpProfil)
            tmpMotif.append(candidatei)
        tmpScore = score(tmpMotif)
        if tmpScore < candidateScore:
            candidateScore = tmpScore
            candidateMotif = tmpMotif
    return candidateMotif, candidateScore


def greedy_motif_search_pseudocount(dnaStr, k, t):
    candidateMotif = [dnaStr[i][0:k] for i in range(len(dnaStr))]
    candidateScore = k * t
    for m0 in range(len(dnaStr[0]) - k + 1):
        tmpMotif = [dnaStr[0][m0:m0 + k]]
        for i in range(1, t):
            tmpProfil = profile_matrix_dict(profile_matrix_pseudocount(tmpMotif))
            candidatei, candidateprobai = get_best_probability_for_profile(dnaStr[i], k, tmpProfil)
            tmpMotif.append(candidatei)
        tmpScore = score(tmpMotif)
        if tmpScore < candidateScore:
            candidateScore = tmpScore
            candidateMotif = tmpMotif
    return candidateMotif, candidateScore


def get_random_profile_from_dna(dnaStr, k):
    #m = len(dnaStr)
    n = len(dnaStr[0])
    result = []
    for i in dnaStr:
        l = random.randint(0, n - k)
        result.append(i[l:l + k])
    return result


def get_best_motif_from_profile(profileDict, k, dnaStr):
    """ finds the best kmer for profile in each strand of dnaStr
    """
    m = len(dnaStr)
    result = []
    for i in range(m):
        tmpBestkmrt, tmpProb = get_best_probability_for_profile(dnaStr[i], k, profileDict)
        result.append(tmpBestkmrt)
    return result


def accumulate(iterable):
    'Return running totals'
    # accumulate([1,2,3,4,5]) --> 1 3 6 10 15
    it = iter(iterable)
    total = next(it)
    yield total
    for element in it:
        total = total + element
        yield total


def biased_random(items, probs):
    accprobs = list(accumulate(probs))
    ubound = accprobs[-1]
    x = random.uniform(0, ubound)
    res = bisect.bisect_right(accprobs, x)

    return items[res]


def randomized_motif_search(dnaStr, k, t):
    """ runs a randomized motif search : select random motif, then tries to improve it  """
    tmpMotif = get_random_profile_from_dna(dnaStr, k)
    bestMotif = tmpMotif
    bestScore = k * t
    end = False
    while not end:
        tmpProfile = profile_matrix_dict(profile_matrix_pseudocount(tmpMotif))
        tmpMotif = get_best_motif_from_profile(tmpProfile, k, dnaStr)
        tmpScore = score(tmpMotif)
        if tmpScore < bestScore:
            bestMotif = tmpMotif
            bestScore = tmpScore
        else:
            end = True
            return bestMotif, bestScore

    exit(-1)


def iterated_randomized_motif_search(dnaStr, k, t, n):
    """ runs randomized_motif_search n times """
    bestMotif = []
    bestScore = k * t

    for i in range(n):
        if i % 10 == 0:
            print "iter :", i
        tmpMotif, tmpScore = randomized_motif_search(dnaStr, k, t)
        if tmpScore < bestScore:
            bestMotif = tmpMotif
            bestScore = tmpScore

    return bestMotif, bestScore


def get_proba_for_all_mker_from_profile(profileDict, k, dnaStrand):
    l = len(dnaStrand)
    resItem = []
    resProba = []
    for i in range(l - k + 1):
        tempItem = dnaStrand[i:i + k]
        tempProba = probability_from_profile(tempItem, profileDict)
        resItem.append(tempItem)
        resProba.append(tempProba)

    return (resItem, resProba)


def gibbs_sampler_search(dnaStr, k, t, n):
    """ runs gibbs sampler """
    tmpMotif = get_random_profile_from_dna(dnaStr, k)
    bestMotif = tmpMotif
    bestScore = score(tmpMotif)

    for m0 in range(n):
        if m0 % 10 == 0:
            print "iter :", m0
        i = random.randint(0, t - 1)
        partialMotif = tmpMotif[:i] + tmpMotif[i + 1:]

        tmpProfile = profile_matrix_dict(profile_matrix_pseudocount(partialMotif))
        (kmlist, probalist) = get_proba_for_all_mker_from_profile(tmpProfile, k, dnaStr[i])
        randomkmer = biased_random(kmlist, probalist)
        tmpMotif[i] = randomkmer
        tmpScore = score(tmpMotif)
        if tmpScore < bestScore:
            bestMotif = tmpMotif
            bestScore = tmpScore

    return bestMotif, bestScore


def iterated_gibbs_sampler_search(dnaStr, k, t, n, rs):
    """ runs randomized_motif_search n times """
    bestMotif = []
    bestScore = k * t

    for i in range(rs):
        print "iter :", i
        tmpMotif, tmpScore = gibbs_sampler_search(dnaStr, k, t, n)
        if tmpScore < bestScore:
            bestMotif = tmpMotif
            bestScore = tmpScore

    return bestMotif, bestScore


def build_prefix_dict(dna_strings):
    res = dict()
    for i in dna_strings:
        tmpprefix = i[:-1]
        if tmpprefix in res:
            res[tmpprefix].append(i)
        else:
            res[tmpprefix] = [i]
    return res


def get_overlap_dict(dna_strings):
    res = dict()
    prefdict = build_prefix_dict(dna_strings)
    for i in dna_strings:
        tmpsuffix = i[1:]
        if tmpsuffix in prefdict:
            res[i] = prefdict[tmpsuffix]
    return res


def build_debruijn_dict(dna_strings):
    res = dict()
    for i in dna_strings:
        tmpprefix = i[:-1]
        tmpsuffix = i[1:]
        if tmpprefix in res:
            res[tmpprefix].append(tmpsuffix)
        else:
            res[tmpprefix] = [tmpsuffix]
    return res


def print_graph(adjacency_list_dict):
    reslist = []
    for i in sorted(sorted(adjacency_list_dict)):
        reslist.append(" ".join([i, "->", ",".join(adjacency_list_dict[i])]))
    res = "\n".join(reslist)
    return res


def find_euler_cycle(adjacency_list_dict):
    adjld = adjacency_list_dict.copy()
    notcomplete = set()
    euler = []
    entrypoint = ''
    exitpoint = ''
    finished = False
    # find a sub graph
    key = adjld.keys()[0]
    euler.append(key)
    while not finished:
        value = adjld[key].pop()
        if len(adjld[key]) == 0:
            adjld.pop(key)
            if key in notcomplete:
                notcomplete.remove(key)
        else:
            notcomplete.add(key)
            #adjld[key] = adjld[key][1:]
        euler.append(value)

        if value in adjld:
            # next step
            key = value
        # we're blocked
        elif len(notcomplete) > 0:
            # take a not complete node to begin new subgraph
            key = notcomplete.pop()
            keypos = euler.index(key)
            euler = euler[keypos:] + euler[1:keypos]
            euler.append(key)
        else:
            # nomore nodes
            if len(adjld) == 0:
                finished = True
            else:
                # Houston, pb
                finished = True
    return euler


def find_euler_path(adjacency_list_dict):
    [startnode, endnode] = find_unbalanced_nodes(adjacency_list_dict)
    if endnode in adjacency_list_dict.keys():
        adjacency_list_dict[endnode].append(startnode)
    else:
        adjacency_list_dict[endnode] = [startnode]

    euler_path = find_euler_cycle(adjacency_list_dict)

    cutpos = 0
    for i in range(len(euler_path) - 1):
        if euler_path[i] == endnode and euler_path[i + 1] == startnode:
            cutpos = i

    euler_path = euler_path[cutpos + 1:] + euler_path[1:cutpos + 1]
    return euler_path


def euler_path_to_nodelist(euler_path):
    res = []
    for i in range(len(euler_path) - 1):
        tmpres = euler_path[i] + euler_path[i + 1][-1]
        res.append(tmpres)
    return res


def adjacency_list_from_string(adj_string):
    adj_list = dict()
    for i in adj_string:
        key, values = i.split(' -> ')
        value = values.split(',')
        adj_list[key] = value
    return adj_list


def out_adjacency_list_from_vmatrix(adj_list, matrix):
    m = len(matrix)
    n = len(matrix[0])
    for i in range(m):
        for j in range(n):
            start = (i, j)
            end = (i + 1, j)
            if end in adj_list:
                adj_list[end].append((start, matrix[i][j]))
            else:
                adj_list[end] = [(start, matrix[i][j])]
    return adj_list


def out_adjacency_list_from_hmatrix(adj_list, matrix):
    m = len(matrix)
    n = len(matrix[0])
    for i in range(m):
        for j in range(n):
            start = (i, j)
            end = (i, j + 1)
            if end in adj_list:
                adj_list[end].append((start, matrix[i][j]))
            else:
                adj_list[end] = [(start, matrix[i][j])]
    return adj_list


def out_adjacency_list_from_diagmatrix(adj_list, matrix):
    m = len(matrix)
    n = len(matrix[0])
    for i in range(m):
        for j in range(n):
            start = (i, j)
            end = (i + 1, j + 1)
            if end in adj_list:
                adj_list[end].append((start, matrix[i][j]))
            else:
                adj_list[end] = [(start, matrix[i][j])]
    return adj_list


def out_adjacency_list_from_strings(strA, strB):
    m = len(strA)
    n = len(strB)
    adj_list = dict()

    for i in range(m, -1, -1):
        for j in range(n, -1, -1):
            curpos = (i, j)
            dpos = (i - 1, j)
            rpos = (i, j - 1)
            diagpos = (i - 1, j - 1)
            if not curpos in adj_list:
                adj_list[curpos] = []
            if i > 0:
                adj_list[curpos].append([dpos, 0])
            if j > 0:
                adj_list[curpos].append([rpos, 0])
            if i > 0 and j > 0:
                if strA[i - 1] == strB[j - 1]:
                    w = 1
                else:
                    w = 0.001
                adj_list[curpos].append([diagpos, w])

    return adj_list


def out_adjacency_list_from_strings_penalty(strA, strB, adj_list):
    m = len(strA)
    n = len(strB)

    for i in range(m, -1, -1):
        for j in range(n, -1, -1):
            curpos = (i, j)
            dpos = (i - 1, j)
            rpos = (i, j - 1)
            diagpos = (i - 1, j - 1)
            if not curpos in adj_list:
                adj_list[curpos] = []
            if j > 0:
                adj_list[curpos].append([rpos, -1])
            if i > 0:
                adj_list[curpos].append([dpos, -1])

            if i > 0 and j > 0:
                if strA[i - 1] == strB[j - 1]:
                    w = 1
                else:
                    w = -1
                adj_list[curpos].append([diagpos, w])

    return adj_list


def out_adjacency_list_from_strings_overlap(strA, strB, adj_list):
    m = len(strA)
    n = len(strB)

    for i in range(m, -1, -1):
        for j in range(n, -1, -1):
            curpos = (i, j)
            dpos = (i - 1, j)
            rpos = (i, j - 1)
            diagpos = (i - 1, j - 1)
            if not curpos in adj_list:
                adj_list[curpos] = []
            if i > 0:
                adj_list[curpos].append([dpos, -2])
            if j > 0:
                adj_list[curpos].append([rpos, -2])
            if i > 0 and j > 0:
                if strA[i - 1] == strB[j - 1]:
                    w = 1
                else:
                    w = -2
                adj_list[curpos].append([diagpos, w])
    return adj_list


def out_adjacency_list_from_strings_with_weight(strA, strB, sigma, wmatrix):
    m = len(strA)
    n = len(strB)
    adj_list = dict()

    for i in range(m, -1, -1):
        for j in range(n, -1, -1):
            curpos = (i, j)
            dpos = (i - 1, j)
            rpos = (i, j - 1)
            diagpos = (i - 1, j - 1)
            if not curpos in adj_list:
                adj_list[curpos] = []

            if j > 0:
                adj_list[curpos].append([rpos, sigma])
            if i > 0:
                adj_list[curpos].append([dpos, sigma])
            if i > 0 and j > 0:
                a1 = strA[i - 1]
                a2 = strB[j - 1]
                w = wmatrix[a1][a2]
                adj_list[curpos].append([diagpos, w])

    return adj_list


def out_fitting_taxi_rides_start(strA, strB, adj_list):
    m = len(strA)
    n = len(strB)

    if not (m, n) in adj_list:
        adj_list[(m, n)] = []
    for i in range(m, 0, -1):
        if not (i, 0) in adj_list:
            adj_list[(i, 0)] = []
        adj_list[(i, 0)].append(((i - 1, 0), 0))

    return adj_list


def out_fitting_taxi_rides_end_fit(strA, strB, adj_list):
    m = len(strA)
    n = len(strB)

    if not (m, n) in adj_list:
        adj_list[(m, n)] = []
    for i in range(m, 0, -1):
        adj_list[(m, n)] = [((i - 1, n), 0)] + adj_list[(m, n)]

    return adj_list


def out_fitting_taxi_rides_start_overlap(strA, strB, adj_list):
    m = len(strA)
    n = len(strB)

    if not (m, n) in adj_list:
        adj_list[(m, n)] = []
    for i in range(m, 0, -1):
        if not (i, 0) in adj_list:
            adj_list[(i, 0)] = []
        adj_list[(i, 0)] = adj_list[(i, 0)] + [((0, 0), 0)]

    for i in range(n, 0, -1):
        if not (0, i) in adj_list:
            adj_list[(0, i)] = []
        adj_list[(0, i)].append(((0, 0), 0))

    return adj_list


def out_fitting_taxi_rides_end_overlap(strA, strB, adj_list):
    m = len(strA)
    n = len(strB)

    if not (m, n) in adj_list:
        adj_list[(m, n)] = []
    for i in range(n, 0, -1):
        adj_list[(m, n)] = adj_list[(m, n)] + [((m, i - 1), 0)]

    for i in range(m, 0, -1):
        adj_list[(m, n)] = adj_list[(m, n)] + [((i - 1, n), 0)]

    return adj_list


def out_taxi_rides(adj_list):
    last_point = max(adj_list.keys())
    for key, value in adj_list.iteritems():
        value.append(((0, 0), 0))
        adj_list[last_point].append((key, 0))

    return adj_list


def print_euler_path(euler_path):
    res = "->".join(euler_path)
    return res


def print_string_from_node_list(node_list):
    res = [node_list[0][:-1]]
    k = len(node_list[0]) - 1
    for i in node_list[:-k]:
        res.append(i[-1])

    #print res
    res = "".join(res)
    return res


def print_string_from_node_list_full(node_list):
    res = [node_list[0][:-1]]
    k = len(node_list[0]) - 1
    for i in node_list:
        res.append(i[-1])

    res = "".join(res)
    return res


def find_unbalanced_nodes(adj_list_dic):
    """ Returns a dict with node as key and a list with count of [in, out]
    """
    inout_counter = dict()
    for key, values in adj_list_dic.iteritems():
        if key in inout_counter:
            inout_counter[key][0] += len(values)
        else:
            inout_counter[key] = [len(values), 0]
        for value in values:
            if value in inout_counter:
                inout_counter[value][1] += 1
            else:
                inout_counter[value] = [0, 1]
    startpoint = ""
    endpoint = ""
    for key, values in inout_counter.iteritems():
        outgoingedges, incomingedges = values
        if outgoingedges > incomingedges:
            startpoint = key
        elif outgoingedges < incomingedges:
            endpoint = key

    return [startpoint, endpoint]


def build_inout_counter(adj_list_dic):
    """ Returns a dict with node as key and a list with count of [in, out]
    """
    inout_counter = dict()
    for key, values in adj_list_dic.iteritems():
        if key in inout_counter:
            inout_counter[key][1] += len(values)
        else:
            inout_counter[key] = [0, len(values)]
        for value in values:
            if value in inout_counter:
                inout_counter[value][0] += 1
            else:
                inout_counter[value] = [1, 0]
    return inout_counter


def get_paired_composition(text, k, d):
    res = []
    for i in range(len(text) - (2 * k) - d + 1):
        kmer1 = text[i:i + k]
        kmer2 = text[i + k + d:i + d + 2 * k]
        res.append([kmer1, kmer2])

    return res


def print_paired_composition(pairedlist):
    tmpres = ["({0}|{1})".format(x[0], x[1]) for x in pairedlist]
    res = ", ".join(sorted(tmpres))

    return res


def print_paired_composition_unsorted(pairedlist):
    tmpres = ["({0}|{1})".format(x[0], x[1]) for x in pairedlist]
    res = ", ".join((tmpres))

    return res


def get_pair_prefix(pairedkmer):
    a = pairedkmer[0][:-1]
    b = pairedkmer[1][:-1]
    return [a, b]


def get_pair_suffix(pairedkmer):
    res1 = pairedkmer[0][1:]
    res2 = pairedkmer[1][1:]
    return [res1, res2]


def get_pair_from_string(string):
    #print string
    return string.split('|')


def build_paired_debruijn_dict(paired_kmer_list):
    l = len(paired_kmer_list)
    res = dict()
    for i in paired_kmer_list:
        tmpprefix = tuple(get_pair_prefix(i))
        tmpsuffix = tuple(get_pair_suffix(i))
        if tmpprefix in res:
            res[tmpprefix].append(tmpsuffix)
        else:
            res[tmpprefix] = [tmpsuffix]
    return res


def print_string_from_paired_node_list(paired_node_list, d):
    node_list = map(lambda x: x[0], paired_node_list)
    node_list_suffix = map(lambda x: x[1], paired_node_list)
    n = len(node_list)
    res = [node_list[0][:-1]]
    k = len(node_list[0])
    #for i in node_list[:-k]:
    for i in node_list:
        res.append(i[-1])
        #res.append('_')
    for i in node_list_suffix[-d - 2:]:
        res.append(i[0])
        #res.append('_')
    res.append(node_list_suffix[-1][1:])
    #print res
    res = "".join(res)
    return res


def build_contigs_from_node(start_node, adj_list_dict, inout_dict):
    contigs_list = []
    for next_node in adj_list_dict[start_node]:
        cur_contig = [start_node]
        loop_node = next_node
        while (inout_dict[loop_node][0] == 1) and (inout_dict[loop_node][1] == 1):
            cur_contig.append(loop_node)
            loop_node = adj_list_dict[loop_node][0]
        cur_contig.append(loop_node)
        contigs_list.append(cur_contig)
    return contigs_list


def get_contigs_from_graph(adj_list_dict):
    inout_dict = build_inout_counter(adj_list_dict)
    contiglist = []
    for starting_node in adj_list_dict.keys():
        inout = inout_dict[starting_node]
        if (inout[0] <> 1) or (inout[1] > 1):
            curcontig = build_contigs_from_node(starting_node, adj_list_dict, inout_dict)
            contiglist += curcontig

    return contiglist


def is_list_compatible_for_turnpike(new_point, point_list, sum_list):
    sum_list_new = sum_list[:]
    for point in (point_list):
        if (point - new_point) in sum_list_new:
            sum_list_new.remove(point - new_point)
            sum_list_new.remove(-point + new_point)
        else:
            return False, []
    sum_list_new.remove(0)
    return True, sum_list_new


def get_sub_turnpike(point_list, sum_list, solution_list):
    if len(sum_list) == 0:
        print "Res", point_list
        return (True, point_list, solution_list)
    new_point_list = point_list[:]
    next_point = max(sum_list)


    # try left
    isok_point, new_sum_list = is_list_compatible_for_turnpike(next_point, point_list, sum_list)
    if isok_point:
        #print "Left" , point_list, next_point
        if len(point_list) == 95:
            print new_sum_list
            print sorted(point_list)
        isok_path, result_path, tmp_sol = get_sub_turnpike(point_list + [next_point], new_sum_list, solution_list)
        if isok_path:
            solution_list += [sorted(result_path)]
            # try right
    next_point = max(point_list) - next_point
    isok_point, new_sum_list = is_list_compatible_for_turnpike(next_point, point_list, sum_list)
    if isok_point:
        #print "Right" ,  point_list, next_point
        #print new_sum_list
        if len(point_list) == 95:
            print new_sum_list
        isok_path, result_path, tmp_sol = get_sub_turnpike(point_list + [next_point], new_sum_list, solution_list)
        if isok_path:
            solution_list += [sorted(result_path)]
    return False, [], solution_list


def recursive_change(money, coins):
    if money == 0:
        return 0
    min_num_coin = sys.maxint
    for coin in coins:
        if money >= coin:
            num_coins = recursive_change(money - coin, coins)
            if (num_coins + 1) < min_num_coin:
                min_num_coin = num_coins + 1
    return min_num_coin


def dynamic_change(money, coins):
    coin_map = [sys.maxint] * (money + 1)
    coin_map[0] = 0
    for value in range(1, money + 1):
        min_num_coin = coin_map[value]
        for coin in coins:
            if coin <= value:
                tmp_num = coin_map[value - coin]
                if tmp_num + 1 < min_num_coin:
                    min_num_coin = tmp_num + 1
        coin_map[value] = min_num_coin

    return coin_map


def build_max_path_len(m, n, adj_list):
    m += 1
    n += 1
    value_map = [[((0, 0), -sys.maxint - 1)] * n for x in range(m)]
    value_map[0][0] = ((0, 0), 0)
    for i in range(m):
        for j in range(n):
            if (i, j) in adj_list:
                for (orig, w) in adj_list[(i, j)]:
                    origx, origy = orig
                    tmpmax = value_map[origx][origy][1] + w
                    if tmpmax > value_map[i][j][1]:
                        value_map[i][j] = (orig, tmpmax)
    return value_map


def build_max_path_len(m, n, adj_list):
    m += 1
    n += 1
    value_map = [[((0, 0), -sys.maxint - 1)] * n for x in range(m)]
    value_map[0][0] = ((0, 0), 0)
    for i in range(m):
        for j in range(n):
            if (i, j) in adj_list:
                for (orig, w) in adj_list[(i, j)]:
                    origx, origy = orig
                    tmpmax = value_map[origx][origy][1] + w
                    if tmpmax > value_map[i][j][1]:
                        value_map[i][j] = (orig, tmpmax)
    return value_map


def print_value_map(value_map):
    values = [[x[1] for x in y] for y in value_map]
    pprint(values)


def backtrack_from_value_map(value_map):
    m = len(value_map)
    n = len(value_map[0])
    curpos = (m - 1, n - 1)
    score = value_map[m - 1][n - 1]
    path = [curpos]

    while curpos != (0, 0):
        posx, posy = curpos
        curpos, score = value_map[posx][posy]
        path.append(curpos)
    path.pop()
    path.reverse()
    path.pop()
    return path


def backtrack_from_value_map2(value_map):
    m = len(value_map)
    n = len(value_map[0])
    endpos = (m - 1, n - 1 )
    posx, posy = endpos
    curpos, score = value_map[posx][posy]

    # if score of last pos == score of previous, then it was a taxi ride !
    posx, posy = curpos
    if score > value_map[posx][posy][1]:
        path = [endpos, curpos]
    else:
        path = [curpos]
    while not ((score == 0) and (curpos == (0, 0))):
        posx, posy = curpos
        curpos, score = value_map[posx][posy]
        path.append(curpos)
    path.pop()
    path.reverse()

    return path


def string_alignement_from_path(strA, strB, path):
    curpos = path[0]
    outstrA = ""
    outstrB = ""
    for nextpos in path[1:]:
        curposx, curposy = curpos
        nextposx, nextposy = nextpos
        if (curposx + 1 == nextposx):
            outstrA += strA[curposx]
        else:
            outstrA += "-"
        if (curposy + 1 == nextposy):
            outstrB += strB[curposy]
        else:
            outstrB += "-"
        curpos = nextpos

    return outstrA, outstrB


def lcs_from_two_aligned_strings(strA, strB):
    lcs = ""
    for i in range(len(strA)):
        if strA[i] == strB[i]:
            lcs += strA[i]
    return lcs


def build_max_path_len_dag(m, adj_list):
    m += 1
    value_map = [(0, -sys.maxint - 1)] * m
    value_map[0] = (0, 0)
    for i in range(m):
        if i in adj_list:
            for (orig, w) in adj_list[i]:
                tmpmax = value_map[orig][1] + w
                if tmpmax > value_map[i][1]:
                    value_map[i] = (orig, tmpmax)
    return value_map


def build_max_path_len_dag_sub(source, sink, adj_list):
    m = sink + 1
    value_map = [(0, -sys.maxint - 1)] * m
    value_map[0] = (0, 0)
    value_map[source] = (source, 0)
    for i in range(source, sink + 1):
        if i in adj_list:
            for (orig, w) in adj_list[i]:
                tmpmax = value_map[orig][1] + w
                if tmpmax > value_map[i][1]:
                    value_map[i] = (orig, tmpmax)
    return value_map


def adjacency_list_from_string_dag(adj_string):
    adj_list = dict()
    adj_list[0] = [(0, 0)]
    for i in adj_string:
        orig, values = i.split('->')
        orig = int(orig)
        dest, weight = values.split(':')
        dest = int(dest)
        weight = int(weight)
        if not dest in adj_list:
            adj_list[dest] = []
        adj_list[dest].append((orig, weight))
    return adj_list


def backtrack_from_value_map_dag(value_map):
    m = len(value_map)
    curpos = (m - 1)
    path = [curpos]
    while curpos != 0:
        curpos = value_map[curpos][0]
        path.append(curpos)
    path.reverse()
    return path


def backtrack_from_value_map_dag_sub(source, sink, value_map):
    curpos = sink
    path = [curpos]
    while curpos != source:
        curpos = value_map[curpos][0]
        path.append(curpos)
        if len(path) > 50:
            print path
            break
    path.reverse()
    return path


def trim_adj_list_dag(source, sink, adj_list):
    for key, value in adj_list.iteritems():
        if key < source or key > sink:
            adj_list.pop(key)
        adj_list[key] = [(i, j) for i, j in value if i >= source and i <= sink]
    return adj_list


def get_blosum_score_windows():
    matrix = {}
    with open("C:\Users\dame\Documents\coursera\w4\BLOSUM62.txt") as f:
        headers = f.readline().split()
        for row in f.readlines():
            elements = row.split()
            matrix[elements[0]] = {headers[index]: int(value) for index, value in enumerate(elements[1:])}
    return matrix


def get_blosum_score():
    matrix = {}
    with open("/home/vador/coursera/w4/BLOSUM62.txt") as f:
        headers = f.readline().split()
        for row in f.readlines():
            elements = row.split()
            matrix[elements[0]] = {headers[index]: int(value) for index, value in enumerate(elements[1:])}
    return matrix


def get_pam_score():
    matrix = {}
    with open("/home/vador/coursera/w4/PAM250_1.txt") as f:
        headers = f.readline().split()
        for row in f.readlines():
            elements = row.split()
            matrix[elements[0]] = {headers[index]: int(value) for index, value in enumerate(elements[1:])}
    return matrix


def find_lev_distance(a, b):
    m = len(a)
    n = len(b)
    m += 1
    n += 1
    value_map = [[( sys.maxint, (0, 0) )] * (n) for x in range(m)]
    value_map[0][0] = (0, (0, 0))
    for i in range(1, m):
        value_map[i][0] = (i, (i - 1, 0))
    for j in range(1, n):
        value_map[0][j] = (j, (0, j - 1))
    for i in range(1, m):
        for j in range(1, n):
            scoreleft = value_map[i - 1][j][0] + 1
            scoreup = value_map[i][j - 1][0] + 1
            scorediag = value_map[i - 1][j - 1][0]
            if a[i - 1] != b[j - 1]:
                scorediag = scorediag + 1
            for valscore, valpoint in [(scorediag, (i - 1, j - 1)), (scoreleft, (i - 1, j)), (scoreup, (i, j - 1))]:
                if valscore < value_map[i][j][0]:
                    value_map[i][j] = (valscore, valpoint)

    return value_map


def find_lev_distance_fit(a, b):
    m = len(a)
    n = len(b)
    m += 1
    n += 1
    value_map = [[( sys.maxint, (0, 0) )] * (n) for x in range(m)]
    value_map[0][0] = (0, (0, 0))
    for i in range(1, m):
        value_map[i][0] = (i, (i - 1, 0))
    for j in range(1, n):
        value_map[0][j] = (j, (0, j - 1))
    for i in range(1, m):
        for j in range(1, n):
            scoreorig = 0
            scoreleft = value_map[i - 1][j][0] + 1
            scoreup = value_map[i][j - 1][0] + 1
            scorediag = value_map[i - 1][j - 1][0]
            if a[i - 1] != b[j - 1]:
                scorediag = scorediag + 1
            for valscore, valpoint in [(scoreorig, (0, 0)), (scorediag, (i - 1, j - 1)), (scoreleft, (i - 1, j)),
                                       (scoreup, (i, j - 1))]:
                if valscore < value_map[i][j][0]:
                    value_map[i][j] = (valscore, valpoint)

    return value_map


def build_max_path_len_opening_penalty(strA, strB, wmatrix):
    sigma = 11
    eps = 1
    m = len(strA)
    n = len(strB)
    value_map_lower = value_map = [[((0, 0, -1), -sys.maxint - 1)] * (n + 1) for x in range(m + 1)]
    value_map_middle = value_map = [[((0, 0, 0), -sys.maxint - 1)] * (n + 1) for x in range(m + 1)]
    value_map_upper = value_map = [[((0, 0, 1), -sys.maxint - 1)] * (n + 1) for x in range(m + 1)]
    value_map = [value_map_lower, value_map_middle, value_map_upper]
    #value_map_middle[0][0] = ((0,0,0),0)
    #value_map_upper[0][0] = ((0,0,0),sigma)
    #value_map_lower[0][0] = ((0,0,0),sigma)

    for i in range(m + 1):
        for j in range(n + 1):

            # fill lower
            # unreachable values if i == 0
            if i > 0:
                (_, scorel) = value_map_lower[i - 1][j]
                value_map_lower[i][j] = ((i - 1, j, -1), scorel - eps)
                (_, scorem) = value_map_middle[i - 1][j]
                if (scorem - sigma) > (scorel - eps):
                    value_map_lower[i][j] = ((i - 1, j, 0), scorem - sigma)

            # Fill upper
            # unreachable values if j == 0
            if j > 0:
                (_, scoreu) = value_map_upper[i][j - 1]
                value_map_upper[i][j] = ((i, j - 1, 1), scoreu - eps)
                (_, scorem) = value_map_middle[i][j - 1]
                if (scorem - sigma) > (scoreu - eps):
                    value_map_upper[i][j] = ((i, j - 1, 0), scorem - sigma)

            # fill middle row
            (_, scorel) = value_map_lower[i][j]
            value_map_middle[i][j] = ((i, j, -1), scorel)

            if i > 0 and j > 0:
                a1 = strA[i - 1]
                a2 = strB[j - 1]
                w = wmatrix[a1][a2]
                (_, scorem) = value_map_middle[i - 1][j - 1]
                scorem += w
                if scorem > scorel:
                    value_map_middle[i][j] = ((i - 1, j - 1, 0), scorem)
                    scorel = scorem

            (_, scoreu) = value_map_upper[i][j]
            if scoreu > scorel:
                value_map_middle[i][j] = ((i, j, 1), scoreu)

            # bootstrap :)
            if (i == 0) and (j == 0):
                value_map_middle[0][0] = ((0, 0, 0), 0)

    return value_map


def backtrack_from_value_map_insertion_penalty(value_map_3):
    m = len(value_map_3[0])
    n = len(value_map_3[0][0])
    value_map = dict()
    value_map[-1] = value_map_3[0]
    value_map[0] = value_map_3[1]
    value_map[1] = value_map_3[2]
    endpos = (m - 1, n - 1, 0 )
    posx, posy, posz = endpos
    curpos, score = value_map[posz][posx][posy]

    # if score of last pos == score of previous, then it was a taxi ride !
    posx, posy, posz = curpos
    if score > value_map[posz][posx][posy][1]:
        path = [endpos, curpos]
    else:
        path = [curpos]
    while not ((score == 0) and (curpos == (0, 0, 0))):
        posx, posy, posz = curpos
        curpos, score = value_map[posz][posx][posy]
        path.append(curpos)
    path.pop()
    path.reverse()

    return path


def string_alignement_from_path_3(strA, strB, path):
    curpos = path[0]
    outstrA = ""
    outstrB = ""
    for nextpos in path[1:]:
        curposx, curposy, curposz = curpos
        nextposx, nextposy, nextposz = nextpos
        if (curposx == nextposx) and (curposy == nextposy):
            pass
        else:
            if (curposx + 1 == nextposx):
                outstrA += strA[curposx]
            else:
                outstrA += "-"
            if (curposy + 1 == nextposy):
                outstrB += strB[curposy]
            else:
                outstrB += "-"
        curpos = nextpos

    return outstrA, outstrB


def linear_path(strA, strB, wmatrix):
    penalty = 5
    m = len(strA)
    n = len(strB)
    curscores = [("-", -sys.maxint - 1) for i in range(m + 1)]
    for j in range(n + 1):
        basescores = curscores
        curscores = [("-", -sys.maxint - 1) for i in range(m + 1)]
        for i in range(m + 1):
            if j == 0:
                if i == 0:
                    curscores[0] = ("-", 0)
                else:
                    _, scoreu = curscores[i - 1]
                    curscores[i] = ("d", scoreu - penalty)
            else:
                _, scorel = basescores[i]
                curscores[i] = ("r", scorel - penalty)
                if i > 0:
                    _, scoreu = curscores[i - 1]
                    if scoreu > scorel:
                        curscores[i] = ("d", scoreu - penalty)
                        scorel = scoreu
                    _, scored = basescores[i - 1]
                    scored += wmatrix[strA[i - 1]][strB[j - 1]]
                    if scored > scorel - penalty:
                        curscores[i] = ("g", scored)
                        scoreu = scorel
                    _, scored = basescores[i - 1]
                    #curscores[i] = max(basescores(i)-penalty,curscores(i-1)-penalty, basescores(i-1)+wmatrix(i,j))
    return curscores


def next_node_from_direction(node, direction):
    posx, posy = node
    if direction == 'd':
        posx += 1
    elif direction == 'r':
        posy += 1
    elif direction == 'g':
        posx += 1
        posy += 1
    return (posx, posy)


def get_middle_egde_linear(strA, strB, matweight):
    n = len(strB)
    ipath1 = linear_path(strA, strB[:n / 2], matweight)
    ipath2 = list(reversed(linear_path(strA[::-1], strB[n / 2:][::-1], matweight)))
    ipaths = [((0, 0), 0) for x in ipath1]
    curval = -sys.maxint - 1
    curidx = -1

    for i in range(len(ipath1)):
        ipaths[i] = (ipath2[i][0], ipath1[i][1] + ipath2[i][1])
        if ipaths[i][1] > curval:
            curval = ipaths[i][1]
            curidx = i

    #print ipaths
    #print  curidx, n/2, ipaths[curidx], curval
    midnode = (curidx, n / 2)
    direction, score = ipaths[curidx]
    nextnode = next_node_from_direction(midnode, direction)
    return midnode, nextnode, direction, score


def linear_space_align(strA, strB, matweight):
    m = len(strA)
    n = len(strB)
    res = ""
    if n == 0:
        #print "d"*m,
        return -5 * m, "d" * m
    midnode, nextnode, direction, score = get_middle_egde_linear(strA, strB, matweight)
    strAB = strA[:midnode[0]]
    strBB = strB[:midnode[1]]
    print(midnode, nextnode)
    print(strAB, "'", strBB)
    _, res = linear_space_align(strAB, strBB, matweight)
    #print midnode, nextnode, direction
    res += direction
    strAC = strA[nextnode[0]:]
    strBC = strB[nextnode[1]:]
    print(strAC, "+", strBC)
    _, tmpres = linear_space_align(strAC, strBC, matweight)
    res += tmpres
    return score, res


def build_str_align_from_directions(strA, strB, directions):
    lstrA = strA[:] + "+++++++"
    lstrB = strB[:] + "+++++++"
    outA = ""
    outB = ""
    ix = 0
    iy = 0
    for d in directions:
        if d == "g":
            outA += lstrA[ix]
            ix += 1
            outB += lstrB[iy]
            iy += 1
        if d == "d":
            outA += lstrA[ix]
            ix += 1
            outB += "-"
        if d == "r":
            outA += '-'
            outB += lstrB[iy]
            iy += 1
    return outA, outB


from multiprocessing import Process, Pipe, Value, Array, Pool


def multiproc_worker_demo(pipe, scalar, array):
    inputval = pipe.recv()
    print('child got : ', inputval)

    n, count, sleeptime = inputval
    scalar.value += 1
    for i in range(len(array)):
        array[i] += 2

    res = [n] * count
    pipe.send(res)
    pipe.close()


def multiproc_master_demo():
    print('coucou')
    (parentEnd, childEnd) = Pipe()
    (parentEnd2, childEnd2) = Pipe()
    print('poy')
    val1 = Value('i', 3)
    val2 = 5
    arr1 = Array('d', 3)
    arr1[0]
    arr2 = [2, 3, 4]
    p1 = Process(target=multiproc_worker_demo, args=(childEnd, val1, arr1))
    p2 = Process(target=multiproc_worker_demo, args=(childEnd2, val2, arr2))
    p1.start()
    p2.start()
    parentEnd.send((10, 20, 5))
    parentEnd2.send((15, 10, 3))
    print('parent got : ', parentEnd.recv())
    print('parent got : ', parentEnd2.recv())
    parentEnd.close()
    parentEnd2.close()
    p1.join()
    p2.join()
    print(val1)
    print(val1)




def linear_space_align_multi(strA, strB, matweight, workerpool):
    m = len(strA)
    n = len(strB)
    res = ""
    if n == 0:
        #print "d"*m,
        return -5 * m, "d" * m
    midnode, nextnode, direction, score = get_middle_egde_linear_multi(strA, strB, matweight, workerpool)
    strAB = strA[:midnode[0]]
    strBB = strB[:midnode[1]]

    _, res = linear_space_align_multi(strAB, strBB, matweight, workerpool)
    #print midnode, nextnode, direction
    res += direction
    strAC = strA[nextnode[0]:]
    strBC = strB[nextnode[1]:]
    _, tmpres = linear_space_align_multi(strAC, strBC, matweight, workerpool)
    res += tmpres
    return score, res


def get_middle_egde_linear_multi(strA, strB, matweight, workerpool):
    n = len(strB)
    args1 = (strA, strB[:n / 2], matweight)
    args2 = (strA[::-1], strB[n / 2:][::-1], matweight)
    (ipath1, ipath2) = workerpool.map(linear_path_multi, [args1, args2])
    ipath2 = list(reversed(ipath2))
    #ipath1 = linear_path(strA, strB[:n / 2], matweight)
    #ipath2 = list(reversed(linear_path(strA[::-1], strB[n / 2:][::-1], matweight)))
    ipaths = [((0, 0), 0) for x in ipath1]
    curval = -sys.maxint - 1
    curidx = -1

    for i in range(len(ipath1)):
        ipaths[i] = (ipath2[i][0], ipath1[i][1] + ipath2[i][1])
        if ipaths[i][1] > curval:
            curval = ipaths[i][1]
            curidx = i

    #print ipaths
    #print  curidx, n/2, ipaths[curidx], curval
    midnode = (curidx, n / 2)
    direction, score = ipaths[curidx]
    nextnode = next_node_from_direction(midnode, direction)
    return midnode, nextnode, direction, score


def get_middle_egde_linear_multi2(strA, strB, matweight):
    n = len(strB)
    ipath1 = linear_path(strA, strB[:n / 2], matweight)
    ipath2 = list(reversed(linear_path(strA[::-1], strB[n / 2:][::-1], matweight)))
    ipaths = [((0, 0), 0) for x in ipath1]
    curval = -sys.maxint - 1
    curidx = -1

    for i in range(len(ipath1)):
        ipaths[i] = (ipath2[i][0], ipath1[i][1] + ipath2[i][1])
        if ipaths[i][1] > curval:
            curval = ipaths[i][1]
            curidx = i

    #print ipaths
    #print  curidx, n/2, ipaths[curidx], curval
    midnode = (curidx, n / 2)
    direction, score = ipaths[curidx]
    nextnode = next_node_from_direction(midnode, direction)
    return midnode, nextnode, direction, score


def linear_path_multi(inparam):
    (strA, strB, wmatrix) = inparam
    penalty = 5
    m = len(strA)
    n = len(strB)
    curscores = [("-", -sys.maxint - 1) for i in range(m + 1)]
    for j in range(n + 1):
        basescores = curscores
        curscores = [("-", -sys.maxint - 1) for i in range(m + 1)]
        for i in range(m + 1):
            if j == 0:
                if i == 0:
                    curscores[0] = ("-", 0)
                else:
                    _, scoreu = curscores[i - 1]
                    curscores[i] = ("d", scoreu - penalty)
            else:
                _, scorel = basescores[i]
                curscores[i] = ("r", scorel - penalty)
                if i > 0:
                    _, scoreu = curscores[i - 1]
                    if scoreu > scorel:
                        curscores[i] = ("d", scoreu - penalty)
                        scorel = scoreu
                    _, scored = basescores[i - 1]
                    scored += wmatrix[strA[i - 1]][strB[j - 1]]
                    if scored > scorel - penalty:
                        curscores[i] = ("g", scored)
                        scoreu = scorel
                    _, scored = basescores[i - 1]
                    #curscores[i] = max(basescores(i)-penalty,curscores(i-1)-penalty, basescores(i-1)+wmatrix(i,j))
    return curscores


def get_3MLCS_path(strA, strB, strC):
    m = len(strA)
    n = len(strB)
    p = len(strC)
    value_map = dict()
    value_map[(0, 0, 0)] = [0, (0, 0, 0)]
    for i in range(m + 1):
        for j in range(n + 1):
            for k in range(p + 1):
                tmpvals = []
                curnode = (i, j, k)
                curscore = -sys.maxint - 1
                prevnode = (0, 0, 0)
                if i > 0:
                    tmpnode = (i - 1, j, k)
                    tmpvals.append((value_map[tmpnode][0], tmpnode))
                if j > 0:
                    tmpnode = (i, j - 1, k)
                    tmpvals.append((value_map[tmpnode][0], tmpnode))
                if k > 0:
                    tmpnode = (i, j, k - 1)
                    tmpvals.append((value_map[tmpnode][0], tmpnode))
                if i > 0 and j > 0:
                    tmpnode = (i - 1, j - 1, k)
                    tmpvals.append((value_map[tmpnode][0], tmpnode))
                if i > 0 and k > 0:
                    tmpnode = (i - 1, j, k - 1)
                    tmpvals.append((value_map[tmpnode][0], tmpnode))
                if j > 0 and k > 0:
                    tmpnode = (i, j - 1, k - 1)
                    tmpvals.append((value_map[tmpnode][0], tmpnode))
                if i > 0 and j > 0 and k > 0:
                    tmpnode = (i - 1, j - 1, k - 1)
                    if strA[i - 1] == strB[j - 1] and strA[i - 1] == strC[k - 1]:
                        w = 1
                    else:
                        w = 0
                    tmpvals.append((value_map[tmpnode][0] + w, tmpnode))
                if i == 0 and j == 0 and k == 0:
                    tmpvals = [[0, (0, 0, 0)]]

                value_map[(i, j, k)] = max(tmpvals, key=(lambda x: x[0]))
    return value_map


def backtrack_from_value_map_3MLCS(value_map):
    sink = max(value_map.keys())
    score = value_map[sink][0]
    curnode = sink
    path = [sink]
    while curnode != (0, 0, 0):
        curnode = value_map[curnode][1]
        path.append(curnode)
    return reversed(path)


def string_alignement_from_path_3MLCS(strA, strB, strC, path):
    curpos = path[0]
    outstrA = ""
    outstrB = ""
    outstrC = ""
    for nextpos in path[1:]:
        curposx, curposy, curposz = curpos
        nextposx, nextposy, nextposz = nextpos
        if (curposx == nextposx) and (curposy == nextposy) and (curposz == nextposz):
            pass
        else:
            if (curposx + 1 == nextposx):
                outstrA += strA[curposx]
            else:
                outstrA += "-"
            if (curposy + 1 == nextposy):
                outstrB += strB[curposy]
            else:
                outstrB += "-"
            if (curposz + 1 == nextposz):
                outstrC += strC[curposz]
            else:
                outstrC += "-"
        curpos = nextpos

    return outstrA, outstrB, outstrC


def print_permutation(permutation):
    print "(" + " ".join(["%+d" % i for i in permutation]) + ")"


def kreversal(source, left, right):
    #print left,right
    if left == right:
        source[left] *= -1
    for i in range((right - left) / 2 + 1):
        #print i, (right-left+1)/2
        tmp = source[left + i]
        source[left + i] = -source[right - i]
        source[right - i] = -tmp
        #source[left+i], source[right-i] = -source[right-i], -source[left+i]
    return None


def greedy_sort_by_reversal(permutation):
    apRevDist = 0
    for i in range(len(permutation)):
        k = i + 1
        if abs(permutation[i]) != k:
            n = permutation.index(k) if k in permutation else permutation.index(-k)
            apRevDist += 1
            kreversal(permutation, i, n)
            print_permutation(permutation)
        if permutation[i] == -k:
            permutation[i] = -permutation[i]
            apRevDist += 1
            print_permutation(permutation)
    return apRevDist


def breakpoints_count(permutation):
    nbBrk = 0
    perm = [0]
    perm.extend(permutation)
    perm.append(len(perm))
    for i in range(len(perm) - 1):
        if (perm[i + 1] - perm[i]) != 1:
            nbBrk += 1
    return nbBrk


def twobreak_count_blocks(strcycles):
    count = 0
    links = dict()
    listcycle = strcycles[1:-1].split(")(")
    for strcycle in listcycle:
        count += len(strcycle.split())
    return count


def twobreak_build_links(strcycles):
    links = dict()
    listcycle = strcycles[1:-1].split(")(")

    for strcycle in listcycle:
        cycle = map(int, strcycle.split())
        origcycle = cycle[0]
        source = origcycle
        for dest in cycle[1:]:
            if source > 0:
                sourcept = (source, 'h')
            else:
                sourcept = (-source, 't')
            if dest > 0:
                destpt = (dest, 't')
            else:
                destpt = (-dest, 'h')
            links[sourcept] = destpt
            links[destpt] = sourcept
            source = dest
        dest = origcycle
        if source > 0:
            sourcept = (source, 'h')
        else:
            sourcept = (-source, 't')
        if dest > 0:
            destpt = (dest, 't')
        else:
            destpt = (-dest, 'h')
        links[sourcept] = destpt
        links[destpt] = sourcept

    return links


def twobreaks_count_cycles(alinks, blinks):
    nblinks = 0
    while len(alinks) > 0:
        nblinks += 1
        curlink = alinks.keys()[0]
        while curlink is not None:
            nextlink = alinks.pop(curlink, None)
            alinks.pop(nextlink, None)
            curlink = blinks.pop(nextlink, None)
            blinks.pop(curlink, None)
    return nblinks


def reversekmer(kmer):
    result = ''
    cc = ''
    for i in range(len(kmer)):
        c = kmer[-i - 1]
        if c == 'A': cc = 'T'
        if c == 'T': cc = 'A'
        if c == 'G': cc = 'C'
        if c == 'C': cc = 'G'
        result += cc
    return result


def built_prefix_dict(dnsStr, k):
    prefix_dict = dict()
    for i in range(len(dnsStr) - k + 1):
        kmer = dnsStr[i:i + k]
        rkmer = reversekmer(kmer)
        if not kmer in prefix_dict:
            prefix_dict[kmer] = []
        prefix_dict[kmer].append(i)
        if rkmer != kmer:
            if not rkmer in prefix_dict:
                prefix_dict[rkmer] = []
            prefix_dict[rkmer].append(i)
    return prefix_dict


def find_shared_kmer(strA, strB, k):
    # put shortest seq in dict
    shared_kmers = []
    invert = False
    if len(strA) < len(strB):
        invert = True
        strA, strB = strB, strA
    p_dict = built_prefix_dict(strB, k)

    for i in range(len(strA) - k + 1):
        kmer = strA[i:i + k]
        if kmer in p_dict:
            for j in p_dict[kmer]:
                shared_kmers.append((i, j))

    if invert:
        shared_kmers = map(lambda (x, y): (y, x), shared_kmers)

    return shared_kmers
