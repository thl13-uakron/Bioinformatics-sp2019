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
    def __init__(self, forward, reverse):
        self.forward = forward
        self.reverse = reverse
        return

    """ fill empty strand with segment of complementary DNA bracketing target subsequence """
    """ assume forward strand of target given """
    def bindPrimer(self, target, length):
        
        if len(self.reverse) == 0: # forward strand present: add reverse primer at back of target
            index = self.forward.find(target) + len(target)
            
            self.reverse = (" " * (index - length)) # empty space
            self.reverse = self.reverse + getComplement(self.forward[index - length : index]) # primer sequence
            self.reverse = self.reverse + (" " * (len(self.forward) - index)) # empty space
            
        elif len(self.forward) == 0: # reverse strand present: add forward primer at front of target
            index = self.reverse.find(getComplement(target))
            
            self.forward = (" " * index) # empty space
            self.forward = self.forward + getComplement(self.reverse[index : index + length]) # primer sequence
            self.forward = self.forward + (" " * (len(self.reverse) - (index + length))) # empty space
            # primers padded with whitespace to keep position identifiable (may change method later)
            
        else: # operation not application with both strands present
            return

    """ extend primer in proper direction to a given length """
    def extendPrimer(self, length):
        # shorter strand identified as primer
        lenDiff = len(self.forward) - len(self.reverse)
        
        if lenDiff > 0: # forward primer extended forward
            # find relative position of primer
            return
        
        elif lenDiff < 0: # reverse primer extended backwards
            # find relative position of primer
            return
        
        else: # operation not applicable in full double-helix
            return

    """put reverse strand into new object, clean up empty spaces, return"""
    def denature(self):
        seperated = DNA("", self.reverse.replace(" ", ""))
        self.forward = self.forward.replace(" ", "")
        self.reverse = ""
        return seperated

    """ display strands """
    def print(self):
        # forward
        print("[5']" + self.forward + "[3']")
        # backward
        print("[5']" + self.reverse + "[3']")


#"""program"""

n = 40 #"""length of original template"""
m = 20 #"""length of segment to amplify"""
p = 5 #"""length of primers"""

d = 200 # base fall-off rate for polymerase
e = 50 # random variation in fall-off rate

cycles = 50

template = generateStrand(n) # original forward dna strand

targetIndex = random.randint(0, n - m - 1) # position of desired dna sequence
forwardTarget = template[targetIndex : targetIndex + m] # desired dna sequence (forward strand)
# reverseTarget = getComplement(forwardTarget) # desired dna sequence (reverse strand)

strands = [DNA(template, getComplement(template))] #"""collection of all dna strands in the simulation"""
strands.append(strands[0].denature())
for s in strands:
    s.bindPrimer(forwardTarget, p)
    s.print()

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

