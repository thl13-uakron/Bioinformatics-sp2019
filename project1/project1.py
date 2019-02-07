# libraries
import random
import re

#"""methods and templates"""

"""create a random string of dna bases of a set length"""
def generateStrand(length):
    bases = ("a", "t", "c", "g")
    DNA = ""
    while len(DNA) < length:
        DNA = DNA + bases[random.randint(0, 3)]
    return DNA

"""check if a DNA string contains only valid base characters"""
def isDNA(string): 
    temp = string.lower()
    temp = temp.replace("a", "")
    temp = temp.replace("g", "")
    temp = temp.replace("c", "")
    temp = temp.replace("t", "")
    return len(temp) == 0

"""get the complementary strand of a DNA string"""
def getComplement(DNA): 
    comp = DNA.lower()
    """match bases, use capital letters as temporary values"""
    comp = comp.replace("a", "T")
    comp = comp.replace("t", "A")
    comp = comp.replace("g", "C")
    comp = comp.replace("c", "G")
    """swap out temporary characters"""
    return comp.lower()

class DNA:
    def __init__():
        #

    """ fill empty strand with segment of complementary DNA bracketing target subsequence """
    def bindPrimer(target, length):
        # forward strand present: add reverse primer at back of target

        # reverse strand present: add forward primer at front of target

        # operation not application with both strands present

    """ extend primer in proper direction to a given length """
    def extendPrimer(length):
        # shorter strand identified as primer

        # forward primer extended forward

        # reverse primer extended backwards

        # operation not applicable in full double-helix

    """ display strands """
    def print(self):
        # forward

        # backward


#"""program"""

n = 2000 #"""length of original template"""
m = 200 #"""length of segment to amplify"""
p = 20 #"""length of primers"""

d = 200 # base fall-off rate for polymerase
e = 50 # random variation in fall-off rate

cycles = 50

template = generateStrand(n) # original forward dna strand

targetIndex = random.randint(0, n - m - 1) # position of desired dna sequence
forwardTarget = template[targetIndex, targetIndex + m] # desired dna sequence (forward strand)
reverseTarget = getComplement(forwardTarget) # desired dna sequence (reverse strand)

strands = [template, getComplement(template)] #"""collection of all dna strands in the simulation"""

"""
PCR Process
- Get subsequence of DNA to replicate
- Perform cycles
    - Match primers to each DNA segment
        - Forward primers match reverse strands, contain first p bases of target
        - Reverse primers match forward strands, contain last p bases of target
    - Extend primers to contain subsequence
        - Length set randomly to d + [-e, e], capped at length of template strand
    - Display intermediate results
    - Repeat with resulting set of strands
- Display final results
"""

