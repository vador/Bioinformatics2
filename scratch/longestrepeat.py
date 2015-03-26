__author__ = 'dame'


import bioalgo.SuffixTree2 as S

a = S.SuffixTree("panamabananas")

class Node():
    def __init__(self, startpos = -1, length = 0):
        self.children = []
        self.startpos = startpos

    def add_child(self, child):
        self.children.append(child)

# Build a BFS
root = []
max_node_id = len(a.nodes)
for i in range(max_node_id):
    root.append(Node())
values = a.edges.values()
values.sort(key=lambda x: x.source_node_index)
for edge in values:


