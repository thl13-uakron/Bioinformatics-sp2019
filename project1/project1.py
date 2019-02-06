"""methods and templates"""

def generateStrand(length): """create a random string of dna bases of a set length"""

def isDNA(string): """check if a DNA string contains only valid base characters"""
    temp = string.lower()
    temp = temp.replace("a", "")
    temp = temp.replace("g", "")
    temp = temp.replace("c", "")
    temp = temp.replace("t", "")
    return len(temp) == 0

def getComplement(DNA): """get the complementary strand of a DNA string"""
    comp = DNA.lower()
    """match bases, use capital letters as temporary values"""
    comp = comp.replace("a", "T")
    comp = comp.replace("t", "A")
    comp = comp.replace("g", "C")
    comp = comp.replace("c", "G")
    """swap out temporary characters"""
    return comp.lower()

def forwardBind(DNA, length): """return a complementary forward primer strand of a given length"""

def reverseBind(DNA, length): """return a complementary reverse primer strand of a given length""""

def forwardExtend(DNA, primer, length): """extend a forward primer to a given length"""

def reverseExtend(DNA, primer, length): """extend a reverse primer to a given length"""


"""program"""

n = 2000 """length of original template"""
m = 200 """length of segment to amplify"""
p = 20 """length of primers"""

strands = [] """collection of all dna strands in the simulation"""
