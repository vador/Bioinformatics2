import math

VERSION = "0.0.1"
__author__ = 'vador'

# A node is a 4 elem list,values being numbers for next node

# ACGT 0123

def convert_nucl_num(strNucl):
    n = len(strNucl)
    res = []
    for char in strNucl:
        if char == 'A':
            res.append(0)
        elif char == 'C':
            res.append(1)
        elif char == 'G':
            res.append(2)
        elif char == 'T':
            res.append(3)
    return res


def convert_num_nucl(listNum):
    res = []
    cvtStr = 'ACGT'
    for i in listNum:
        res.append(cvtStr[i])
    return "".join(res)


def get_next_trie_pos(val, sourcenode):
    return sourcenode[val]


def add_sub_to_trie(trie, lstValues):
    repeat = 0
    curnodepos = 0
    for value in lstValues:
        curnode = trie[curnodepos]
        # if already in trie :
        if curnode[value] is not None:
            curnodepos = curnode[value]
            repeat += 1
        else:
            # must create the node !
            trie.append([None, None, None, None])
            curnodepos = len(trie) - 1
            curnode[value] = curnodepos
    return repeat


def dump_trie_edges(trie):
    cvtStr = 'ACGT'
    for i in range(len(trie)):
        node = trie[i]
        for j in range(4):
            next = node[j]
            if next is not None:
                print i + 1, next + 1, cvtStr[j]


def is_node_leaf(node):
    return not any(node)


def prefix_trie_matching(num_text, trie):
    i = 0
    s = num_text[i]
    curnode = trie[0]
    pattern = []
    while True:
        if curnode[s] is not None and i < len(num_text):
            curnode = trie[curnode[s]]
            pattern.append(s)
            i += 1
            if i < len(num_text):
                s = num_text[i]
        elif is_node_leaf(curnode):
            return pattern
        else:
            return None


def trie_matching(num_text, trie):
    for i in range(len(num_text) - 1):
        if prefix_trie_matching(num_text[i:], trie) is not None:
            yield (i)


def suffix_array(strTxt):
    return sf_arr_qsort([x for x in range(len(strTxt))], strTxt)


def partial_suffix_array(strTxt, k):
    suf_arr = suffix_array(strTxt)
    res = []
    for i in range(len(suf_arr)):
        if suf_arr[i] % k == 0:
            res.append((i, suf_arr[i]))
    return res


def cmp_suffixes(strTxt, i1, i2):
    # strings are equal !
    io1 = i1
    io2 = i2
    if i1 == i2:
        return 0
        # else, they must be different so it will end before out of range!
    n = len(strTxt)
    s1 = strTxt[i1]
    s2 = strTxt[i2]
    while True:
        if s1 < s2:
            return 1
        elif s1 > s2:
            return -1
        else:
            i1 += 1
            i2 += 1
            if i2 >= len(strTxt):
                print strTxt
                print io1, io2, i2
                print strTxt[io1], strTxt[io2]
            s1 = strTxt[i1]
            s2 = strTxt[i2]


def sf_arr_qsort(list, strTxt):
    """Quicksort using list comprehensions"""
    if list == []:
        return []
    else:
        pivot = list[0]
        lorg = [cmp_suffixes(strTxt, pivot, x) for x in list[1:]]
        tmplist1 = [list[x] for x in range(1, len(list)) if lorg[x - 1] < 0]
        tmplist2 = [list[x] for x in range(1, len(list)) if lorg[x - 1] >= 0]
        lesser = sf_arr_qsort(tmplist1, strTxt)
        greater = sf_arr_qsort(tmplist2, strTxt)
        return lesser + [pivot] + greater

