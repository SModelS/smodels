"""
.. module:: tools.summaryReader
   :synopsis: Classes to read the summary.txt files.
    
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    

"""

class Output():
    def __init__(self, l):
        self.ana = l[0]
        self.topo = l[1]
        self.sqrts = l[2]
        self.condViolation = l[3]
        self.tval = l[4]
        self.ul = l[5]
        self.rval = l[6]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False


class Summary():
    """
    Class to access the output given in the summary.txt
    
    """
    def __init__(self, filename):
        self.results = []
        self.read(filename)

    def read(self, infile):
        f = open(infile)
        resultLines = None
        for l in f:
            if not l.strip():
                continue
            if "#Analysis" in l:
                resultLines = True
                continue
            if l.startswith("#"):
                continue
            if "----" in l:
                continue
            if "====" in l:
                resultLines = None
            if resultLines: self.results.append(Output(l.split()))

    def __eq__(self, other):
        if not len(self.results) == len(other.results):
            return False
        for res in self.results:
            if not res in other.results:
                return False
        return True
