""" rna2d module holds functions to use for RNA 2D representation, comparison,
and some generic operations"""


class RNA2dStructure(object):
    """2D structure representation of RNA
    stores only information about base pairs in Vienna format"""
    __slots__ = ['dot_bracket', 'energy', 'identification', 'seq']

    def __init__(self,
                 dot_bracket,
                 energy=float("inf"),
                 identification="",
                 seq=""):
        assert is_valid_dot_bracket(dot_bracket)
        self.dot_bracket = dot_bracket
        self.energy = energy
        self.identification = identification
        self.seq = check_format_rna(seq)

    def __str__(self):
        return self.dot_bracket + "\n" + self.seq

    def __repr__(self):
        return self.__str__()


def only_paired(dot_bracket):
    """remove '.' character from dot_bracket"""
    assert is_valid_dot_bracket(dot_bracket), "invalid dot-bracket"
    return "".join([i for i in dot_bracket if i != '.'])


def check_format_rna(seq):
    """verifies the identity of the nucleotides in the RNA sequence"""
    result = [char.upper() for char in seq]
    for nucleotide in result:
        assert nucleotide in ['A', 'T', 'C', 'G', 'U'], "illegal character"
        if nucleotide == 'T':
            nucleotide = 'U'
    return result


def is_valid_dot_bracket(dot_bracket):
    """tests Vienna dot-bracket for illegal structure (or symbol)"""
    counter = 0
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


def dot_bracket_to_mountain(dot_bracket):
    """Vienna dotbracket to mountain representation
    e.g. "..(((.....))).." -> [0, 0, 1, 2, 3, 3, 3, 3, 3, 3, 2, 1, 0, 0, 0]"""
    assert is_valid_dot_bracket(dot_bracket), "invalid dot-bracket"
    counter = 0
    val = list()
    for i in dot_bracket:
        if(i == '('):
            counter += 1
        elif(i == ')'):
            counter -= 1
        val.append(counter)
    return val


def dot_bracket_to_bp_set(dot_bracket):
    """Vienna dotbracket to base pair set
    (sorted list by first base of the pair)
    "((..))" -> [(1,6), (2,5)]"""
    bpset = list()
    accumulator = list()
    count = 1  # 1 indexed
    for i in dot_bracket:
        if (i == '('):
            accumulator.append(count)
        elif(i == ')'):
            bpset.append((accumulator.pop(), count))
        count += 1
    bpset.sort()
    return bpset


def mountain_distance(mountain1, mountain2):
    """lp1 mountain distance on two mountains representation of same length
    e.g. [1,2,2,2,1], [1,2,3,2,1] = 1"""
    assert len(mountain1) == len(mountain2)

    def absdiff(num_tuple):
        """absolute difference, applied over a zipped array"""
        return abs(num_tuple[0] - num_tuple[1])

    return sum([absdiff(num_tuple) for num_tuple in zip(mountain1, mountain2)])


def base_pair_set_distance(base_pair1, base_pair2):
    """base pair set distance  |(A - B) U (B - A)|"""
    return len((set(base_pair1).symmetric_difference(set(base_pair2))))


def to_bpseq(rna_sequence, dot_bracket):
    """outputs RNA sequence and structure to bpseq format
    index base paired (0->unpaired)
    1 C 0
    2 C 9
    3 U 8
    4 G 0
    5 A 0
    6 A 0
    7 C 0
    8 A 3
    9 G 2"""
    assert is_valid_dot_bracket(dot_bracket), "invalid dot-bracket"
    assert len(rna_sequence) == len(dot_bracket)
    rna_sequence = check_format_rna(rna_sequence)
    base_pairs = dot_bracket_to_bp_set(dot_bracket)
    result = list()  # (Int, Char, Int)
    result.append("header")
    for (i, base_pair) in enumerate(rna_sequence):
        result.append((i+1, base_pair, 0))
    for (base_pair1, base_pair2) in base_pairs:
        result[base_pair1] = (result[base_pair1][0],
                              result[base_pair1][1],
                              base_pair2)
        result[base_pair2] = (result[base_pair2][0],
                              result[base_pair2][1],
                              base_pair1)
    return "\n".join([str(x[0]) + " " + str(x[1]) +
                      " " + str(x[2]) for x in result[1:]])


def get_stems(dot_bracket):
    """gets the stems from the dot_bracket
    ([opening], [closing], dict of pairs"""
    assert is_valid_dot_bracket(dot_bracket), "invalid dot-bracket"
    list_opener = []
    stems = []
    list_stem_end = []
    i = 0
    # separate into stems
    while(i <= len(dot_bracket)-1):
        # add to opening list
        if(dot_bracket[i] == '('):
            list_opener.append(i)
        # find the closing in the opening list
        elif(dot_bracket[i] == ')'):
            stem = ([], [], dict())
            while(i <= len(dot_bracket)-1):
                if(dot_bracket[i] == ')'):
                    opener = list_opener.pop()
                    stem[0].append(opener)
                    stem[1].append(i)
                    stem[2][opener] = i
                    if((not len(list_opener) == 0) and
                       (list_opener[-1] in list_stem_end)):
                        list_stem_end.append(list_opener[-1])
                        break
                elif(dot_bracket[i] == '('):
                    if (not len(list_opener) == 0):
                        list_stem_end.append(list_opener[-1])
                    i -= 1
                    break
                i += 1
            stems.append(stem)
        i += 1
        stems.sort(reverse=True)
    return stems


__all__ = ["RNA2dStructure",
           "check_format_rna",
           "is_valid_dot_bracket",
           "to_bpseq",
           "get_stems"]
