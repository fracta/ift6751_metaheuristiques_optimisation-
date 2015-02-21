"""Node class, used for tree structure"""

import numpy as np
cimport numpy as np

import collections


cdef inline bint richcmp_helper(int compare, int op):
    """Returns True/False for each compare operation given an op code.
    Compare should act similarly to Java's comparable interface"""
    if op == 2: # ==
        return compare == 0
    elif op == 3: # !=
        return compare != 0
    elif op == 0: # <
        return compare < 0
    elif op == 1: # <=
        return compare <= 0
    elif op == 4: # >
        return compare > 0
    elif op == 5: # >=
        return compare >= 0


cdef class Node:
    """ A simple tree representation"""
    cdef public Node parent
    cdef public str label
    cdef public list children
    cdef public int _id

    def __init__(self, Node parent, str label="o", int _id = -1):
        if parent is not None:
            assert isinstance(parent, Node)
            parent.children.append(self)
        self.label = label
        self.children = list()
        self.parent = parent
        self._id = _id

    cpdef list get_children(self):
        """get the children list of the node"""
        return self.children

    cpdef str get_label(self):
        """get label of the node"""
        return self.label

    def __richcmp__(Node self, Node other, int op):
      """used to compare individuals for the ranking in the hall of fame"""
      cdef double v1 = self._id
      cdef double v2 = other._id
      if v1 > v2:
          compare = 1
      elif v1 < v2:
          compare = -1
      else:
          compare = 0
      return richcmp_helper(compare, op)

    def  __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.pretty_print_id()

    def pretty_print_id(self, prefix="", is_tail=True):
        """see http://stackoverflow.com/questions/4965335"""
        cdef str result = prefix + (
            "└── " if is_tail else "├── ") + str(self._id) + "\n"

        for index in range(0, len(self.children)-1):
            result += self.children[index].pretty_print_id(
                prefix + ("    " if is_tail else "│   "), False)
        if len(self.children) >= 1:
            result += self.children[-1].pretty_print_id(
                prefix + ("    " if is_tail else "│   "), True)
        return result

    cpdef str pretty_print(self, prefix="", is_tail=True):
        """see http://stackoverflow.com/questions/4965335"""
        cdef str result = prefix + (
            "└── " if is_tail else "├── ") + self.label + "\n"

        for index in range(0, len(self.children)-1):
            result += self.children[index].pretty_print(
                prefix + ("    " if is_tail else "│   "), False)
        if len(self.children) >= 1:
            result += self.children[-1].pretty_print(
                prefix + ("    " if is_tail else "│   "), True)
        return result

    cpdef append(self, Node other_node):
        """add the other_node and the end of the node's children list"""
        self.children.append(other_node)


cpdef bint is_valid_dot_bracket(dot_bracket):
    """tests Vienna dot-bracket for illegal structure (or symbol)"""
    cdef int counter = 0

    for i in dot_bracket:
        if i == '(':
            counter += 1
        elif i == ')':
            counter -= 1
        elif i != '.':  # illegal symbol
            return False
        if counter < 0:  # unbalanced structure
            return False
    if counter != 0:
        return False  # unbalanced structure
    return True


cpdef Node dot_bracket_to_tree(dot_bracket):
    """creates a abstract shape base pair tree from the Vienna dot bracket"""
    assert is_valid_dot_bracket(dot_bracket), "invalid dot-bracket"
    root = Node(None, label='r')
    position = root
    for c in dot_bracket:
        if c == '(':
            child = Node(position, label='l')
            position = child
        elif c == ')':
            position = position.parent
        else:
            continue
    return root


###############################################################################
# PREORDER AND POSTORDER TRAVERSAL LABELING

cdef preorder_helper(Node tree, int* index_ptr):
    """recursive preorder traversal"""
    # perform action before
    tree._id = index_ptr[0]
    index_ptr[0] += 1
    for children in tree.get_children():
        preorder_helper(children, index_ptr)
    return

cpdef label_preorder(Node root):
    cdef int index = 0
    cdef int* index_ptr = &index
    preorder_helper(root, index_ptr)
    return


cdef postorder_helper(Node tree, int* index_ptr):
    """recursive postorder traversal"""
    for children in tree.get_children():
        postorder_helper(children, index_ptr)
    # perform action after
    tree._id = index_ptr[0]
    index_ptr[0] += 1  # update the index value
    return


cpdef label_postorder(Node root):
    cdef int index = 0
    cdef int* index_ptr = &index
    postorder_helper(root, index_ptr)
    return


cpdef list get_ancestors(Node node):
    cdef list result = []
    while node.parent != None:
        result.append(node.parent._id)
        node = node.parent  # move up
    return result

cpdef int get_lmd(Node node):
    """get the id of the leftmost descendent"""
    while len(node.get_children())!= 0:
        node = node.get_children()[0]
    return node._id

cdef postorder_list_helper(Node node, list res):
    """helper for the postorder_list function"""
    for children in node.get_children():
        postorder_list_helper(children, res)
    res.append(node)
    return

cpdef list postorder_list(Node root):
    """returns an array of the nodes in postorder traversal order"""
    cdef list result = list()
    postorder_list_helper(root, result)
    return result



cdef class Tree:
    """ tree object used for the computation of tree edit distance"""
    cdef readonly Node root
    cdef public np.ndarray nodes
    cdef public np.ndarray lmds
    cdef public np.ndarray keyroots

    def __init__(self, root):
        """build the tree with annotation for zhang-shasha algorithm"""

        # root of the tree
        self.root = root

        # postorder emueration of the nodes
        self.nodes = np.array(postorder_list(root))

        # left most descendents
        self.lmds = np.zeros(len(self.nodes), dtype=int)

        # postorder label the tree
        label_postorder(root)

        # figure out the leftmost descendents for each nodes
        for node in self.nodes:
            self.lmds[node._id] = get_lmd(node)

        # figure out the keyroots, {(max k | lmd(k) = x) forall x}
        tmp = collections.defaultdict(lambda:-1)
        for node in reversed(self.nodes):
            index = self.lmds[node._id]
            if tmp[index] == -1:
                tmp[index] = node._id
        self.keyroots = np.array(sorted(tmp.values()))
        return


    def __getitem__(self, index):
        assert index in range(0, len(self.nodes)), "index out of range"
        return self.nodes[index]



###############################################################################
# ZHANG SHASHA
###############################################################################

def simple_distance(A, B, get_children=Node.get_children,
                    get_label=Node.get_label, label_dist=strdist):
    """Computes the exact tree edit distance between trees A and B.
    Use this function if both of these things are true:
    * The cost to insert a node is equivalent to ``label_dist('', new_label)``
    * The cost to remove a node is equivalent to ``label_dist(new_label, '')``
    Otherwise, use :py:func:`zss.distance` instead.
    :param A: The root of a tree.
    :param B: The root of a tree.
    :param get_children:
        A function ``get_children(node) == [node children]``.  Defaults to
        :py:func:`zss.Node.get_children`.
    :param get_label:
        A function ``get_label(node) == 'node label'``.All labels are assumed
        to be strings at this time. Defaults to :py:func:`zss.Node.get_label`.
    :param label_distance:
        A function
        ``label_distance((get_label(node1), get_label(node2)) >= 0``.
        This function should take the output of ``get_label(node)`` and return
        an integer greater or equal to 0 representing how many edits to
        transform the label of ``node1`` into the label of ``node2``. By
        default, this is string edit distance (if available). 0 indicates that
        the labels are the same. A number N represent it takes N changes to
        transform one label into the other.
    :return: An integer distance [0, inf+)
    """
    return distance(
        A, B, get_children,
        insert_cost=lambda node: label_dist('', get_label(node)),
        remove_cost=lambda node: label_dist(get_label(node), ''),
        update_cost=lambda a, b: label_dist(get_label(a), get_label(b)),
    )


#def unlabeled_distance(A, B, get_children=Node.get_children):
    #""" unlabeled tree distance """
    #INS = lambda x: 1
    #DEL = lambda x: 1
    #SUBS = lambda x, y: 0
    #return distance(A, B, get_children, INS, DEL, SUBS)


#def distance(A, B, get_children, insert_cost, remove_cost, update_cost):
    #'''Computes the exact tree edit distance between trees A and B with a
    #richer API than :py:func:`zss.simple_distance`.
    #Use this function if either of these things are true:
    #* The cost to insert a node is **not** equivalent to the cost of changing
      #an empty node to have the new node's label
    #* The cost to remove a node is **not** equivalent to the cost of changing
      #it to a node with an empty label
    #Otherwise, use :py:func:`zss.simple_distance`.
    #:param A: The root of a tree.
    #:param B: The root of a tree.
    #:param get_children:
        #A function ``get_children(node) == [node children]``.  Defaults to
        #:py:func:`zss.Node.get_children`.
    #:param insert_cost:
        #A function ``insert_cost(node) == cost to insert node >= 0``.
    #:param remove_cost:
        #A function ``remove_cost(node) == cost to remove node >= 0``.
    #:param update_cost:
        #A function ``update_cost(a, b) == cost to change a into b >= 0``.
    #:return: An integer distance [0, inf+)
    #'''
    #A, B = Tree(A, get_children), Tree(B, get_children)
    #treedists = np.zeros((len(A.nodes), len(B.nodes)), dtype=int)

    #def treedist(i, j):
        #Al = A.lmds
        #Bl = B.lmds
        #An = A.nodes
        #Bn = B.nodes

        #m = i - Al[i] + 2
        #n = j - Bl[j] + 2
        #fd = np.zeros((m, n), dtype=int)

        #ioff = Al[i] - 1
        #joff = Bl[j] - 1

        #for x in range(1, m):  # δ(l(i1)..i, θ) = δ(l(1i)..1-1, θ) + γ(v → λ)
            #fd[x][0] = fd[x-1][0] + remove_cost(An[x-1])
        #for y in range(1, n):  # δ(θ, l(j1)..j) = δ(θ, l(j1)..j-1) + γ(λ → w)
            #fd[0][y] = fd[0][y-1] + insert_cost(Bn[y-1])

        #for x in range(1, m):  # the plus one is for the xrange impl
            #for y in range(1, n):
                ## only need to check if x is an ancestorsestor of i
                ## and y is an ancestorsestor of j
                #if Al[i] == Al[x+ioff] and Bl[j] == Bl[y+joff]:
                    ##                   +-
                    ##                   | δ(l(i1)..i-1, l(j1)..j) + γ(v → λ)
                    ## δ(F1 , F2 ) = min-+ δ(l(i1)..i , l(j1)..j-1) + γ(λ → w)
                    ##                   | δ(l(i1)..i-1, l(j1)..j-1) + γ(v → w)
                    ##                   +-
                    #fd[x][y] = min(
                        #fd[x-1][y] + remove_cost(An[x+ioff]),
                        #fd[x][y-1] + insert_cost(Bn[y+joff]),
                        #fd[x-1][y-1] + update_cost(An[x+ioff], Bn[y+joff]),
                    #)
                    #treedists[x+ioff][y+joff] = fd[x][y]
                #else:
                    ##                   +-
                    ##                   | δ(l(i1)..i-1, l(j1)..j) + γ(v → λ)
                    ## δ(F1 , F2 ) = min-+ δ(l(i1)..i , l(j1)..j-1) + γ(λ → w)
                    ##                   | δ(l(i1)..l(i)-1, l(j1)..l(j)-1)
                    ##                   |                     + treedist(i1,j1)
                    ##                   +-
                    #p = Al[x+ioff]-1-ioff
                    #q = Bl[y+joff]-1-joff
                    #fd[x][y] = min(
                        #fd[x-1][y] + remove_cost(An[x+ioff]),
                        #fd[x][y-1] + insert_cost(Bn[y+joff]),
                        #fd[p][q] + treedists[x+ioff][y+joff]
                    #)
    #for i in A.keyroots:
        #for j in B.keyroots:
            #treedist(i, j)
    #return treedists[-1][-1]