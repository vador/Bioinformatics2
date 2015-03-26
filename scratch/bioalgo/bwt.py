VERSION = "0.0.1"
__author__ = 'vador'

from collections import deque


def cmp_bwt(strTxt, i1, i2):
    # strings are equal !
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
            i1 = i1 % n
            i2 = i2 % n
            s1 = strTxt[i1]
            s2 = strTxt[i2]


def bwt_arr_qsort(list, strTxt):
    """Quicksort using list comprehensions"""
    if list == []:
        return []
    else:
        pivot = list[0]
        lorg = [cmp_bwt(strTxt, pivot, x) for x in list[1:]]
        tmplist1 = [list[x] for x in range(1, len(list)) if lorg[x - 1] < 0]
        tmplist2 = [list[x] for x in range(1, len(list)) if lorg[x - 1] >= 0]
        lesser = bwt_arr_qsort(tmplist1, strTxt)
        greater = bwt_arr_qsort(tmplist2, strTxt)
        return lesser + [pivot] + greater


def bwt_transform(strTxt):
    bwt_list = bwt_arr_qsort([x for x in range(len(strTxt))], strTxt)
    bwt_res = [strTxt[i - 1] for i in bwt_list]
    return "".join(bwt_res)


def make_base_counter(bwtTxt):
    alphabet = set(bwtTxt)
    base = dict()
    for letter in alphabet:
        base[letter] = 0
    counter_list = [base]
    return counter_list


def make_counter(bwtTxt):
    n = len(bwtTxt)
    counter_list = make_base_counter(bwtTxt)
    for letter in bwtTxt:
        newcount = counter_list[-1].copy()
        newcount[letter] += 1
        counter_list.append(newcount)
    return counter_list


def bwt_add_letter(partial_list, first_col, counter_list):
    init_counter = dict()


def rle_enc(strTxt):
    res = []
    cur_char = ""
    cur_count = 0
    for char in strTxt:
        if char != cur_char:
            res.append((cur_char, cur_count))
            cur_char = char
            cur_count = 0
        cur_count += 1
    res.append((cur_char, cur_count))
    return res[1:]


def reverse_bwt_transform(bwtStr):
    firstcol = sorted(bwtStr)
    lastcol = bwtStr[:]
    n = len(bwtStr)

    for i in range(n):
        tmpcol = []
        #for x in range(n):
        #    tmpcol.append("".join([lastcol[x], firstcol[x]]))
        #firstcol = sorted(tmpcol)
        #firstcol = [lastcol[x]+ firstcol[x] for x in range(n)]
        firstcol = [lc + fc for lc, fc in zip(lastcol, firstcol)]
        firstcol.sort()

        if (i % 100) == 0:
            print i

    result = firstcol[0][1:]
    return result


def build_alphabet_counter(strTxt):
    cnt = dict()
    alphabet = set(strTxt)
    for l in alphabet:
        cnt[l] = []
    for i in range(len(strTxt)):
        cnt[strTxt[i]].append(i)
    return cnt


def build_first_col_cnt(firstcol):
    cnt = dict()
    alphabet = set(firstcol)
    firstcol = "".join(firstcol)
    for l in alphabet:
        cnt[l] = firstcol.find(l)
    return cnt


def reverse_bwt_transform2(bwtStr):
    n = len(bwtStr)
    firstcol = sorted(bwtStr)
    lastcol = bwtStr[:]
    fcnt = build_first_col_cnt(firstcol)
    lcnt = build_alphabet_counter(lastcol)
    resStr = ""
    letter = "$"
    pos = 0

    for i in range(n):
        idx = lcnt[letter][pos]
        letter = firstcol[idx]
        resStr += letter
        pos = idx - fcnt[letter]

    return resStr


def last_to_first(firstcol, lastcol):
    l2f = []
    lcnt = dict()
    alphabet = set(firstcol)
    fcnt = build_first_col_cnt(firstcol)
    for i in alphabet:
        lcnt[i] = 0

    for i in range(len(firstcol)):
        letter = lastcol[i]
        pos = lcnt[letter]
        l2f.append(fcnt[letter] + pos)
        lcnt[letter] += 1

    return l2f


def bwt_matching(firstcol, lastcol, pattern, last2firstcol):
    top = 0
    n = len(lastcol)
    bottom = len(lastcol) - 1
    while top <= bottom:
        if len(pattern) > 0:
            sym = pattern[-1]
            pattern = pattern[:-1]
            #tmpstr = lastcol[top:bottom]
            topindex = n - 1
            bottomindex = 0
            bFound = False
            tmpstr = lastcol[top:bottom + 1]
            for i in range(top, bottom + 1):
                if lastcol[i] == sym:
                    bFound = True
                    if i < topindex:
                        topindex = i
                    bottomindex = i
            if bFound:
                top = last2firstcol[topindex]
                bottom = last2firstcol[bottomindex]
            else:
                return 0
        else:
            return bottom - top + 1


def bwt_better_matching(firstcol, lastcol, pattern, fcnt, lcnt):
    top = 0
    n = len(lastcol)
    bottom = len(lastcol) - 1
    while top <= bottom:
        if len(pattern) > 0:
            sym = pattern[-1]
            pattern = pattern[:-1]

            if lcnt[top][sym] != lcnt[bottom + 1][sym]:
                top = fcnt[sym] + lcnt[top][sym]
                bottom = fcnt[sym] + lcnt[bottom + 1][sym] - 1
            else:
                return 0
        else:
            return bottom - top + 1


def bwt_better_matching_with_pos(firstcol, lastcol, pattern, fcnt, lcnt):
    top = 0
    n = len(lastcol)
    bottom = len(lastcol) - 1
    while top <= bottom:
        if len(pattern) > 0:
            sym = pattern[-1]
            pattern = pattern[:-1]

            if lcnt[top][sym] != lcnt[bottom + 1][sym]:
                top = fcnt[sym] + lcnt[top][sym]
                bottom = fcnt[sym] + lcnt[bottom + 1][sym] - 1
            else:
                return 0
        else:
            return (top, bottom)