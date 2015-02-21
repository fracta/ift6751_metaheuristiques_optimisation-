from tree import *
from rna2d import *

def dot_bracket_to_tree(dot_bracket):
    """creates a abstract shape base pair tree from the Vienna dot bracket"""
    assert is_valid_dot_bracket(dot_bracket), "invalid dot-bracket"
    root = Node(None)
    position = root
    for c in dot_bracket:
        if c == '(':
            child = Node(position)
            position.append(child)
            position = child
        elif c == ')':
            position = position.parent
        else:
            continue
    return root