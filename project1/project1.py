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

"""return a complementary forward primer strand of a given length"""
def forwardBind(DNA, length):
    return

"""return a complementary reverse primer strand of a given length"""
def reverseBind(DNA, length):
    return

"""extend a forward primer to a given length"""
def forwardExtend(DNA, primer, length):
    return

"""extend a reverse primer to a given length"""
def reverseExtend(DNA, primer, length):
    return


#"""program"""

n = 2000 #"""length of original template"""
m = 200 #"""length of segment to amplify"""
p = 20 #"""length of primers"""

template = generateStrand(n) # original dna strand 
strands = [template, getComplement(template)] #"""collection of all dna strands in the simulation"""

