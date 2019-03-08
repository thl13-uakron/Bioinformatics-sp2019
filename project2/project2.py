"""Data Structures"""

"""Helper Functions"""

# Smith-Waterman algorithm
# Given a query and target sequence of nucleotide bases,
# and a scoring matrix:
# (1) Create score grid corresponding to strand bases
# (2) Fill squares ([0, m], 0) and (0, [0, n])
# (3) Fill remaining squares based on edge values
#   (a) Moving down: match gap with s2[j]
#   (b) Moving across: match s1[i] with gap
#   (c) Moving diagonally: match s1[i] with s2[j]
# (4) Get path of optimal score at square (n, m)
# (5) Translate path into alignment, return aligned pair

def fillGrid(s1, s2, scoringMatrix, alignmentGrid):
    # to fill rest of grid:
    # (1) fill square [k][k] starting with k=1
    # (2) extend values in same column
    # (3) extend values in same row
    # (4) incremement k
    # (5) repeat until k reaches length of shorter string                                                                      
    
    # to fill individual squares:
    # (1) find highest score among available options:
    #   (a) moving down from square [i][j - 1]
    #   (b) moving right from square [i - 1][j]
    #   (c) moving diagonally from square [i - 1][j - 1]

    #
    
    return alignmentGrid

def smithWaterman(s1, s2, scoringMatrix):
    # get initial grid values
    alginmentGrid = swInitializeGrid(s1, s2, scoringMatrix)

    # fill remaining values and find optimal path
    alignmentGrid = fillGrid(s1, s2, scoringMatrix, alignmentGrid)
    return swGetTraceback(s1, s2, scoringMatrix, alignmentGrid)

def swInitializeGrid(s1, s2, scoringMatrix):
    # fill initial grid value
    alignmentGrid = [[0]]

    # fill leading row
    for i in range len(s1):
        alignmentGrid[i + 1] = [alignmentGrid[i] + scoringMatrix[s1[i]][" "]]

    # fill leading column
    for j in range len(s2):
        alignmentGrid[0][j + 1] = alignmentGrid[0][j] + scoringMatrix[" "][s2[j]

    # return initial grid
    return alignmentGrid

def swGetTraceback(s1, s2, scoringMatrix, alignmentGrid):
    # To find global alignment pattern from completed grid:
    # (1) Start at square [n][m] at bottom right
    # (2) Compare score with squares [n - 1][m], [n][m - 1], and [n - 1][m - 1]
    # (3) Find option where score difference matches score from corresponding pairing
    # (4) Add bases to aligned pair
    # (5) Repeat from new square until reaching [0][0]
    # (6) Return aligned pair
                                                                           
    result = ("", "")
    
    return result

# BLASTN algorithm
# Given a query and target sequence of nucleotide bases,
# a scoring matrix, query length, and score cutoff:
# (1) Identify and filter out low-complexity regions
# (2) Create a list of sufficiently-scoring sub-queries
# (3) Search for matches of each query in the target sequence
# (4)

"""Driver Program"""

bases = ["A", "T", "C", "G", " "]
scoringMatrix = {"A":{"A":1, "T":0, "G":0, "C":0, " ":-1},
                 "T":{"T":1, "G":0, "C":0, " ":-1},
                 "C":{"C":1, "G":0, " ":-1},
                 "G":{"G":1, " ":-1},
                 " ":{" ":0}}


