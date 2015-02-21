# -*- coding: utf-8 -*-
"""tree edit distancestorse and related structures"""


import numpy as np
cimport numpy as np

import collections
from tree import Tree



def strdist(a, b):
    """simple string equality distancestorse"""
    if a == b:
        return 0
    else:
        return 1


cdef class AnnotatedTree:
    """ tree object used for the computation of tree edit distancestorse"""
    cdef Tree root
    cdef list nodes
    cdef list lmds
    cdef list keyroots

    def __init__(self, root):
        self.root = root
        self.nodes = list()  # a pre-order enumeration of the nodes in the tree
        self.lmds = list()   # left most descendents
        self.keyroots = None
        # k and k' are nodes specified in the pre-order enumeration.
        # keyroots = {k | there exists no k'>k such that lmd(k) == lmd(k')}
        # see paper for more on keyroots
        cdef list stack = list()
        cdef list pstack = list()

        # start the annotation loop
        cdef int j = 0
        stack.append((root, collections.deque()))

        # preorder label the nodes
        while len(stack) > 0:
            node, ancestors = stack.pop()
            print node
            node._id = j
            for children in node.get_children():
                ancestor = collections.deque(ancestors)
                ancestor.appendleft(node._id)
                stack.append((children, ancestor))
            pstack.append((node, ancestors))
            j += 1

        # figure out which are the leftmost descendents
        lmds = dict()
        keyroots = dict()
        cdef int i = 0
        while len(pstack) > 0:
            node, ancestors = pstack.pop()
            self.nodes.append(node)
            if len(node.get_children()) == 0:
                lmd = i
                for ancestor in ancestors:
                    if ancestor not in lmds:
                        lmds[ancestor] = i
                    else:
                        break
            else:
                try:
                    lmd = lmds[node._id]
                except:
                    import pdb
                    pdb.set_trace()
            self.lmds.append(lmd)
            keyroots[lmd] = i
            i += 1
        self.keyroots = sorted(keyroots.values())


def simple_distancestorse(A, B, get_children=Tree.get_children,
                    get_label=Tree.get_label, label_dist=strdist):
    """Computes the exact tree edit distancestorse between trees A and B.
    Use this function if both of these things are true:
    * The cost to insert a node is equivalent to ``label_dist('', new_label)``
    * The cost to remove a node is equivalent to ``label_dist(new_label, '')``
    Otherwise, use :py:func:`zss.distancestorse` instead.
    :param A: The root of a tree.
    :param B: The root of a tree.
    :param get_children:
        A function ``get_children(node) == [node children]``.  Defaults to
        :py:func:`zss.Tree.get_children`.
    :param get_label:
        A function ``get_label(node) == 'node label'``.All labels are assumed
        to be strings at this time. Defaults to :py:func:`zss.Tree.get_label`.
    :param label_distancestorse:
        A function
        ``label_distancestorse((get_label(node1), get_label(node2)) >= 0``.
        This function should take the output of ``get_label(node)`` and return
        an integer greater or equal to 0 representing how many edits to
        transform the label of ``node1`` into the label of ``node2``. By
        default, this is string edit distancestorse (if available). 0 indicates that
        the labels are the same. A number N represent it takes N changes to
        transform one label into the other.
    :return: An integer distancestorse [0, inf+)
    """
    return distancestorse(
        A, B, get_children,
        insert_cost=lambda node: label_dist('', get_label(node)),
        remove_cost=lambda node: label_dist(get_label(node), ''),
        update_cost=lambda a, b: label_dist(get_label(a), get_label(b)),
    )


def unlabeled_distancestorse(A, B, get_children=Tree.get_children):
    """ unlabeled tree distancestorse """
    INS = lambda x: 1
    DEL = lambda x: 1
    SUBS = lambda x, y: 0
    return distancestorse(A, B, get_children, INS, DEL, SUBS)


def distancestorse(A, B, get_children, insert_cost, remove_cost, update_cost):
    '''Computes the exact tree edit distancestorse between trees A and B with a
    richer API than :py:func:`zss.simple_distancestorse`.
    Use this function if either of these things are true:
    * The cost to insert a node is **not** equivalent to the cost of changing
      an empty node to have the new node's label
    * The cost to remove a node is **not** equivalent to the cost of changing
      it to a node with an empty label
    Otherwise, use :py:func:`zss.simple_distancestorse`.
    :param A: The root of a tree.
    :param B: The root of a tree.
    :param get_children:
        A function ``get_children(node) == [node children]``.  Defaults to
        :py:func:`zss.Tree.get_children`.
    :param insert_cost:
        A function ``insert_cost(node) == cost to insert node >= 0``.
    :param remove_cost:
        A function ``remove_cost(node) == cost to remove node >= 0``.
    :param update_cost:
        A function ``update_cost(a, b) == cost to change a into b >= 0``.
    :return: An integer distancestorse [0, inf+)
    '''
    A, B = AnnotatedTree(A, get_children), AnnotatedTree(B, get_children)
    treedists = np.zeros((len(A.nodes), len(B.nodes)), dtype=int)

    def treedist(i, j):
        Al = A.lmds
        Bl = B.lmds
        An = A.nodes
        Bn = B.nodes

        m = i - Al[i] + 2
        n = j - Bl[j] + 2
        fd = np.zeros((m, n), dtype=int)

        ioff = Al[i] - 1
        joff = Bl[j] - 1

        for x in range(1, m):  # δ(l(i1)..i, θ) = δ(l(1i)..1-1, θ) + γ(v → λ)
            fd[x][0] = fd[x-1][0] + remove_cost(An[x-1])
        for y in range(1, n):  # δ(θ, l(j1)..j) = δ(θ, l(j1)..j-1) + γ(λ → w)
            fd[0][y] = fd[0][y-1] + insert_cost(Bn[y-1])

        for x in range(1, m):  # the plus one is for the xrange impl
            for y in range(1, n):
                # only need to check if x is an ancestorsestor of i
                # and y is an ancestorsestor of j
                if Al[i] == Al[x+ioff] and Bl[j] == Bl[y+joff]:
                    #                   +-
                    #                   | δ(l(i1)..i-1, l(j1)..j) + γ(v → λ)
                    # δ(F1 , F2 ) = min-+ δ(l(i1)..i , l(j1)..j-1) + γ(λ → w)
                    #                   | δ(l(i1)..i-1, l(j1)..j-1) + γ(v → w)
                    #                   +-
                    fd[x][y] = min(
                        fd[x-1][y] + remove_cost(An[x+ioff]),
                        fd[x][y-1] + insert_cost(Bn[y+joff]),
                        fd[x-1][y-1] + update_cost(An[x+ioff], Bn[y+joff]),
                    )
                    treedists[x+ioff][y+joff] = fd[x][y]
                else:
                    #                   +-
                    #                   | δ(l(i1)..i-1, l(j1)..j) + γ(v → λ)
                    # δ(F1 , F2 ) = min-+ δ(l(i1)..i , l(j1)..j-1) + γ(λ → w)
                    #                   | δ(l(i1)..l(i)-1, l(j1)..l(j)-1)
                    #                   |                     + treedist(i1,j1)
                    #                   +-
                    p = Al[x+ioff]-1-ioff
                    q = Bl[y+joff]-1-joff
                    fd[x][y] = min(
                        fd[x-1][y] + remove_cost(An[x+ioff]),
                        fd[x][y-1] + insert_cost(Bn[y+joff]),
                        fd[p][q] + treedists[x+ioff][y+joff]
                    )
    for i in A.keyroots:
        for j in B.keyroots:
            treedist(i, j)
    return treedists[-1][-1]


__all__ = ['unlabeled_distancestorse', 'Tree']