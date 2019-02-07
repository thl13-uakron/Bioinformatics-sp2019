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

"""abstraction of DNA pair with PCR operations as methods"""
class DNA:
    def __init__(self, forward, reverse):
        self.forward = {"bases": forward, "pos": 0}
        self.reverse = {"bases": reverse, "pos": 0}
        return

    """ fill empty strand with segment of complementary DNA bracketing target subsequence """
    """ assume forward strand of target given """
    def bindPrimer(self, forwardPrimer, reversePrimer):
        
        if len(self.reverse["bases"]) == 0: # forward strand present: add reverse primer at back of target
            index = self.forward["bases"].find(getComplement(reversePrimer)) # search for primer complement

            if index > -1: # add primer if complement found
                self.reverse["pos"] = index # log primer position
                self.reverse["bases"] = reversePrimer # log primer sequence
            
        elif len(self.forward["bases"]) == 0: # reverse strand present: add forward primer at front of target
            index = self.reverse["bases"].find(getComplement(forwardPrimer)) # search for primer complement

            if index > -1: # add primer if complement found
                self.forward["pos"] = index # log primer position
                self.forward["bases"] = forwardPrimer # log primer sequence
            
        else: # operation not applicable with both strands present
            return

    """ extend primer in proper direction to a given length """
    def extendPrimer(self, length):
        # shorter strand identified as primer
        forwardLen = len(self.forward["bases"])
        reverseLen = len(self.reverse["bases"])
        lenDiff = forwardLen - reverseLen

        # operation not applicable if primer is empty
        if forwardLen == 0 and reverseLen == 0:
            return
        
        if lenDiff < 0: # forward primer extended forward
            primer = self.forward
            template = self.reverse
            
            # cap length as needed to match template
            if template["pos"] + length > reverseLen:
                length = reverseLen - template["pos"]
                
            # extend primer to length towards end of template
            primer["bases"] = getComplement(template["bases"][primer["pos"] : primer["pos"] + length])
            
        
        elif lenDiff > 0: # reverse primer extended backwards
            primer = self.reverse
            template = self.forward
            
            # cap length as needed to match template
            if length > primer["pos"] + len(primer["bases"]):
                length = primer["pos"] + len(primer["bases"])
                
            # extend primer to length towards start of template
            primer["pos"] -= length - len(primer["bases"])
            primer["bases"] = getComplement(template["bases"][primer["pos"] : primer["pos"] + length])
        
        else: # operation not applicable in full double-helix
            return

    """put reverse strand into new object, reset relative positioning, return"""
    def denature(self):
        seperated = DNA("", self.reverse["bases"])
        self.forward["pos"] = 0
        self.reverse = {"bases": "", "pos": 0}
        return seperated

    """ display strands """
    def print(self):
        # forward
        print((" " * self.forward["pos"]) + "5' " + self.forward["bases"] + " 3'")
        # backward
        print((" " * self.reverse["pos"]) + "3' " + self.reverse["bases"] + " 5'")


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

"""split collection of DNA pairs into collection individual strands"""
def denatureAll(strandList):
    newList = []
    for s in strandList: # get seperated strands
        # skip seperation for single strands
        if len(s.forward["bases"]) > 0 and len(s.reverse["bases"]) > 0:
            newList.append(s.denature())
            
    # add new strands back into container
    for n in newList:
        strandList.append(n)

"""add complementary primers to all single strands at target subsequence for replication"""
def bindAll(strandList, forwardTarget, primerLen, variation=0):
    forwardPrimer = forwardTarget[: primerLen]
    reversePrimer = getComplement(forwardTarget[0 - primerLen :])
    
    for s in strandList:
        s.bindPrimer(forwardPrimer, reversePrimer)

"""extend all DNA primers to length in variation range"""
def extendAll(strandList, extensionLen, variation=0):
    for s in strandList:
        s.extendPrimer(extensionLen + random.randint(0 - variation, variation))

"""display contents of all DNA pairs"""
def printAll(strandList):
    print("\n")
    for s in strandList:
        print("")
        s.print()
        
    
#"""driver program"""

n = 80 #"""length of original template"""
m = 20 #"""length of segment to amplify"""
p = 10 #"""length of primers"""

d = 25 # base fall-off rate for polymerase
e = 5 # random variation in fall-off rate

cycles = 10

template = generateStrand(n) # original forward dna strand

targetIndex = random.randint(0, n - m - 1) # position of desired dna sequence
forwardTarget = template[targetIndex : targetIndex + m] # desired dna sequence (forward strand)
# reverseTarget = getComplement(forwardTarget) # desired dna sequence (reverse strand)

strands = [DNA(template, getComplement(template))] #"""collection of all dna strands in the simulation"""

# perform cycles
printAll(strands)
while cycles > 0:
    denatureAll(strands)
    bindAll(strands, forwardTarget, p)
    extendAll(strands, d, e)
    printAll(strands)
    cycles -= 1

