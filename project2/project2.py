"""Data Structures"""

# character to represent gap values
gap = "-"

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

# full algorithm combining all steps
def smithWaterman(s1, s2, scoringMatrix):
    # get initial grid values
    alginmentGrid = swInitializeGrid(s1, s2, scoringMatrix)

    # fill remaining values and find optimal path
    alignmentGrid = swFillGrid(s1, s2, scoringMatrix, alignmentGrid)
    return swGetTraceback(s1, s2, scoringMatrix, alignmentGrid)

# components of algorithm
def swFillGrid(s1, s2, scoringMatrix, alignmentGrid):
    # to fill rest of grid:
    # (1) fill square [k][k] starting with k=1
    # (2) extend values in same column
    # (3) extend values in same row
    # (4) incremement k
    # (5) repeat until k reaches length of shorter string m                                                                     
    
    if len(s1) < len(s2):
        m = len(s1) + 1
    else:
        m = len(s2) + 1

    for k in range (1, m):
        swFillSquare(s1, s2, scoringMatrix, alignmentGrid, k, k) # (1)
        for i in range (k, len(s1)):
            swFillSquare(s1, s2, scoringMatrix, alignmentGrid, i + 1, k) # (2)
        for j in range (k, len(s2)):
            swFillSquare(s1, s2, scoringMatrix, alignmentGrid, k, j + 1) # (3)
        # (4) (5)
    
    return alignmentGrid

def swFillSquare(s1, s2, scoringMatrix, alignmentGrid, x, y):
    # to fill individual squares:
    # (1) find highest score among available options:
    #   (a) moving down from square [i][j - 1]
    #   (b) moving right from square [i - 1][j]
    #   (c) moving diagonally from square [i - 1][j - 1]

    alignmentGrid[x][y] = alignmentGrid[x][y - 1] + scoringMatrix[gap][s2[y - 1]] # (1a)

    score = alignmentGrid[x - 1][y] + scoringMatrix[s1[x - 1]][gap] # (1b)
    if score > alignmentGrid[x][y]:
        alignmentGrid[x][y] = score

    score = alignmentGrid[x - 1][y - 1] + scoringMatrix[s1[x - 1]][s2[y - 1]] # (1c)
    if score > alignmentGrid[x][y]:
        alignmentGrid[x][y] = score

    return 

def swInitializeGrid(s1, s2, scoringMatrix):
    # fill initial grid value
    alignmentGrid = [[0]]

    # fill leading row, create empty columns
    for i in range (0, len(s1)):
        alignmentGrid.append([alignmentGrid[i][0] + scoringMatrix[s1[i]][gap]])
        for j in range (0, len(s2)):
            alignmentGrid[i + 1].append(None)

    # fill leading column
    for j in range (0, len(s2)):
        alignmentGrid[0].append(alignmentGrid[0][j] + scoringMatrix[gap][s2[j]])

    # return initial grid
    return alignmentGrid

def swGetTraceback(s1, s2, scoringMatrix, alignmentGrid):
    # To find global alignment pattern from completed grid:
    # (1) Start at square [n][m] at bottom right
    # (2) Compare score with squares [n - 1][m], [n][m - 1], and [n - 1][m - 1]
    # (3) Find option where score difference matches score from corresponding pairing
    # (4) Add bases to aligned pair, move to new square
    # (5) Repeat from new square until reaching [0][0]
    # (6) Return aligned pair
                                                                           
    result = {s1:"", s2:""}

    # (1)
    n = len(s1)
    m = len(s2)

    while n > 0 or m > 0:
        # (2)
        if n > 0 and alignmentGrid[n][m] - alignmentGrid[n - 1][m] == scoringMatrix[s1[n - 1]][gap]: # (3)
            # (4)
            n -= 1
            result[s1] = s1[n] + result[s1]
            result[s2] = gap + result[s2]
            
        elif m > 0 and alignmentGrid[n][m] - alignmentGrid[n][m - 1] == scoringMatrix[gap][s2[m - 1]]:
            m -= 1
            result[s1] = gap + result[s1]
            result[s2] = s2[m] + result[s2]
            
        elif n > 0 and m > 0 and alignmentGrid[n][m] - alignmentGrid[n - 1][m - 1] == scoringMatrix[s1[n - 1]][s2[m - 1]]:
            n -= 1
            m -= 1
            result[s1] = s1[n] + result[s1]
            result[s2] = s2[m] + result[s2]
        # (5)
    
    return result # (6)

def printGrid(s1, s2, alignmentGrid):
    padding = 3

    row = ("[" + (" " * padding) + "] ") * 2
    for i in range (0, len(s1)):
        row += "[" + (" " * (padding - 1)) + s1[i] + "] "
    print(row)
    
    for j in range (0, len(alignmentGrid[0])):
        row = ""

        if j > 0:
            columnChar = s2[j - 1]
        else:
            columnChar = " "
        row += "[" + (" " * (padding - 1)) + columnChar + "] "
        
        for i in range (0, len(alignmentGrid)):
            val = str(alignmentGrid[i][j])
            row += "[" + (" " *  (padding - len(val))) + val + "] "
        print(row)
        
    return

def swAnalyseAlignment(results):
    # Given returned output from swGetTraceback method
    return

# BLASTN algorithm
# Given a query and target sequence of nucleotide bases,
# a scoring matrix, query length, and score cutoff:
# (1) Identify and filter out low-complexity regions
# (2) Create a list of sufficiently-scoring sub-queries
# (3) Search for matches of each query in the target sequence
# (4)

"""Driver Program"""

# score system
scoringMatrix = {"A":{"A":1, "T":0, "G":0, "C":0, gap:-1},
                 "T":{"T":1, "G":0, "C":0, gap:-1},
                 "C":{"C":1, "G":0, gap:-1},
                 "G":{"G":1, gap:-1},
                 gap:{gap:0}}
for b1 in scoringMatrix:
     for b2 in scoringMatrix[b1]:
         scoringMatrix[b2][b1] = scoringMatrix[b1][b2]

# strand sequences
s1 = "AAACGGGTTTCCCCCAAA"
s2 = "AGCCCCCAAAGTCGTCCCCGTCCGTC"

grid = swInitializeGrid(s1, s2, scoringMatrix)
swFillGrid(s1, s2, scoringMatrix, grid)
printGrid(s1, s2, grid)
traceback = swGetTraceback(s1, s2, scoringMatrix, grid)
print(traceback[s1])
print(traceback[s2])

