import string

import bioalgo.bioalgorithm as B

POSITIVE_INFINITY = 1 << 30


class Node:
    def __init__(self, suffix_link=None):
        self.suffix_link = suffix_link

    def __repr__(self):
        return 'Node(' + str(self.suffix_link) + ')'


class Edge:
    def __init__(self, src_node_idx, dst_node_idx, first_char_idx, last_char_idx):
        self.src_node_idx = src_node_idx
        self.dst_node_idx = dst_node_idx
        self.first_char_idx = first_char_idx
        self.last_char_idx = last_char_idx

    def split(self, suffix, suffix_tree):
        return split_edge(self, suffix, suffix_tree)

    def __len__(self):
        return self.last_char_idx - self.first_char_idx + 1

    def __repr__(self):
        return 'Edge(' + str(self.src_node_idx) + ', ' + \
               str(self.dst_node_idx) + ', ' + \
               str(self.first_char_idx) + ', ' + \
               str(self.last_char_idx) + ')'


def split_edge(edge, suffix, suffix_tree):
    #alloc new node
    new_node = Node()#suffix.src_node_idx
    suffix_tree.nodes.append(new_node)
    new_node_idx = len(suffix_tree.nodes) - 1
    #alloc new edge
    new_edge = Edge(new_node_idx, edge.dst_node_idx, edge.first_char_idx + len(suffix), edge.last_char_idx)
    suffix_tree.insert_edge(new_edge)
    #shorten existing edge
    edge.last_char_idx = edge.first_char_idx + len(suffix) - 1
    edge.dst_node_idx = new_node_idx
    return new_node_idx


class Suffix:
    def __init__(self, src_node_idx, first_char_idx, last_char_idx):
        self.src_node_idx = src_node_idx
        self.first_char_idx = first_char_idx
        self.last_char_idx = last_char_idx

    def is_explicit(self):
        return is_explicit_suffix(self)

    def is_implicit(self):
        return is_implicit_suffix(self)

    def canonize(self, suffix_tree):
        canonize_suffix(self, suffix_tree)

    def __repr__(self):
        return 'Suffix(' + str(self.src_node_idx) + ', ' + str(self.first_char_idx) + ', ' + str(
            self.last_char_idx) + ')'

    def __len__(self):
        return self.last_char_idx - self.first_char_idx + 1


def is_explicit_suffix(suffix):
    return suffix.first_char_idx > suffix.last_char_idx


def is_implicit_suffix(suffix):
    return not is_explicit_suffix(suffix)


def canonize_suffix(suffix, suffix_tree):
    if not suffix.is_explicit():
        edge = suffix_tree.edge_lookup[suffix.src_node_idx, suffix_tree.string[suffix.first_char_idx]]
        if (len(edge) <= len(suffix)):
            suffix.first_char_idx += len(edge)
            suffix.src_node_idx = edge.dst_node_idx
            canonize_suffix(suffix, suffix_tree)


class SuffixTree:
    def __init__(self, string, alphabet=None):
        self.string = string
        if alphabet == None:
            alphabet = set(string)
        self.alphabet = alphabet
        self.nodes = [Node()]
        self.edge_lookup = {} #edge_source_node_first_char_dict
        self.active_point = Suffix(0, 0, -1)
        for i in range(len(string)):
            add_prefix(i, self.active_point, self)

    def insert_edge(self, edge):
        self.edge_lookup[edge.src_node_idx, self.string[edge.first_char_idx]] = edge

    def remove_edge(self, edge):
        del self.edge_lookup[edge.src_node_idx, self.string[edge.first_char_idx]]


def add_prefix(last_char_idx, active_point, suffix_tree):
    last_parent_node_idx = -1
    while True:
        parent_node_idx = active_point.src_node_idx
        if active_point.is_explicit():
            if (
            active_point.src_node_idx, suffix_tree.string[last_char_idx]) in suffix_tree.edge_lookup: #already in tree
                break
        else:
            edge = suffix_tree.edge_lookup[active_point.src_node_idx, suffix_tree.string[active_point.first_char_idx]]
            if suffix_tree.string[edge.first_char_idx + len(active_point)] == suffix_tree.string[
                last_char_idx]: #the given prefix is already in the tree, do nothing
                break
            else:
                parent_node_idx = edge.split(active_point, suffix_tree)
        suffix_tree.nodes.append(Node(-1))
        new_edge = Edge(parent_node_idx, len(suffix_tree.nodes) - 1, last_char_idx, POSITIVE_INFINITY)##################
        suffix_tree.insert_edge(new_edge)
        #add suffix link
        if last_parent_node_idx > 0:
            suffix_tree.nodes[last_parent_node_idx].suffix_link = parent_node_idx
        last_parent_node_idx = parent_node_idx
        if active_point.src_node_idx == 0:
            active_point.first_char_idx += 1
        else:
            active_point.src_node_idx = suffix_tree.nodes[active_point.src_node_idx].suffix_link
        active_point.canonize(suffix_tree)
    if last_parent_node_idx > 0:
        suffix_tree.nodes[last_parent_node_idx].suffix_link = parent_node_idx
        #last_parent_node_idx = parent_node_idx
    active_point.last_char_idx += 1
    active_point.canonize(suffix_tree)


def show_edge(suffix_tree, src_node_idx, first_char):
    edge = suffix_tree.edge_lookup[src_node_idx, first_char]
    #print edge
    print suffix_tree.string[edge.first_char_idx:edge.last_char_idx + 1]


def common_length(suffix_tree, strin):
    cnodeidx = 0
    cl = 0
    curnode = suffix_tree.nodes[0]
    s = strin[cl]
    while cl < len(strin):
        try:
            s = strin[cl]
            curedge = suffix_tree.edge_lookup[cnodeidx, s]
            sbase = suffix_tree.string[curedge.first_char_idx:curedge.last_char_idx + 1]
            if sbase == strin[cl:cl + len(sbase)]:
                cnodeidx = curedge.dst_node_idx
                cl += curedge.last_char_idx - curedge.first_char_idx + 1
            else:
                i = 0
                while (i < len(sbase)) and (cl + i < len(strin)) and (sbase[i] == strin[cl + i]):
                    i += 1
                cl += i
                break
        except KeyError:
            break

    return cl

