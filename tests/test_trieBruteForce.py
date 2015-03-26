from unittest import TestCase
from nose2 import *

from TrieConstruction import TrieBruteForce
__author__ = 'dame'


class TestTrieBruteForce(TestCase):

    def test_basicTrie(self):
        exp = """0->1:A
1->2:T
2->3:A
3->4:G
4->5:A
2->6:C
0->7:G
7->8:A
8->9:T
"""
        myTrie = TrieBruteForce()
        myTrie.addPattern("ATAGA")
        myTrie.addPattern("ATC")
        myTrie.addPattern("GAT")
        res = myTrie.dumpTrie()
        self.assertEquals(exp, res)


    def test_prefixMatch(self):
        text = "AATCGGGTTCAATCGGGGTA"
        myTrie = TrieBruteForce()
        myTrie.addPattern("ATCG")
        myTrie.addPattern("GGGT")

        for i in range(len(text)):
            if (myTrie.prefixMatch(text[i:])):
                print i
