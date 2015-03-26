#!/usr/local/bin/python
# coding: utf-8

__author__ = 'dame'

"""
TRIECONSTRUCTION(Patterns)
        Trie ← a graph consisting of a single node root
        for each string Pattern in Patterns
            currentNode ← root
            for i ← 1 to |Pattern|
                if there is an outgoing edge from currentNode with label currentSymbol
                    currentNode ← ending node of this edge
                else
                    add a new node newNode to Trie
                    add a new edge from currentNode to newNode with label currentSymbol
                    currentNode ← newNode
        return Trie
"""


class TrieBruteForce:
    """ """
    root = None
    class TrieNode:
        label = None
        children = None
        def __init__(self, label = None):
            self.label = label
            self.children = []

        def addChild(self, child):
            self.children.append(child)

        def hasChildren(self):
            return (len(self.children)>0)

        def getChildWithLabel(self, label):
            for child in self.children:
                if child.label == label:
                    return child
            return None

    def __init__(self):
        self.root = self.TrieNode()


    def getChildWithLabel(self, node, label):
        if node is None:
            return None
        for child in node.children:
            if child.label == label:
                return child
        return None


    def addMissingSuffix(self, suffix):
        None ;

    def addPattern(self, pattern, startNode = root):
        if len(pattern) < 1:
            return
        if startNode is None:
            startNode = self.root

        tmpLabel = pattern[0]
        child = self.getChildWithLabel(startNode, tmpLabel)
        if child is None:
            child = self.TrieNode()
            child.label = tmpLabel
            startNode.addChild(child)
        self.addPattern(pattern[1:], child)

    def dumpSubTrie(self,startNode, parentcnt, maxcnt, strres):
        curcnt = parentcnt + 1
        for child in startNode.children:
            strres = strres + str(parentcnt) + "->" + str(curcnt) + ":" + child.label + "\n"
            (maxcnt,strres) = self.dumpSubTrie(child,curcnt, maxcnt,strres)
            curcnt = maxcnt

        return (curcnt, strres)


    def dumpTrie(self):
        strres = ""
        (n, strres) = self.dumpSubTrie(self.root, 0,0, strres)
        return strres

    def prefixMatch(self, text):
        i = 0
        symbol = text[i]
        v = self.root
        pattern = ""
        while True:
            if not v.hasChildren():
                return pattern
            else:
                if i == len(text):
                    return
                symbol = text[i]
                w = v.getChildWithLabel(symbol)
                if w is not None:
                    v = w
                    pattern += symbol
                    i +=1
                else:
                    return
    def trieMatching(self, text):
        for i in range(len(text)):
            if (myTrie.prefixMatch(text[i:])):
                print i

if __name__ == '__main__':
    import sys
    myTrie = TrieBruteForce()
    with open(sys.argv[1]) as f:
        content = f.readlines()
    text = content[0]
    for string in content[1:]:
        string = string.rstrip('\n\r')
        if len(string)>0:
            myTrie.addPattern(string)
    #output = myTrie.dumpTrie()
    myTrie.trieMatching(text)



