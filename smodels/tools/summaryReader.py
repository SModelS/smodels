"""
.. module:: tools.summaryReader
   :synopsis: Classes to read the summary.txt files.

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""

class Output():
    def __init__(self, l, signalRegion, txnames):
        self.ana = l[0]
        self.sqrts = l[1]
        self.condViolation = l[2]
        self.tval = l[3]
        self.ul = l[4]
        self.rval = l[5]
        self.txnames = txnames
        self.signalRegion = signalRegion

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

class MissOutput():
    def __init__(self, l):
        self.sqrts = eval(l[0])
        self.weight = eval(l[1])
        self.element = l[3]

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
        self.signalRegions = []
        self.txnames = []
        self.filename = filename
        self.missedTopos = []
        self.read(filename)

    def __str__ ( self ):
        import os
        fn = self.filename.replace("//","/")
        final = fn.replace ( os.getcwd(), "." )
        if len(final)>50:
            final="..."+final[-47:]
        return "Summary(%s)" % final

    def read(self, infile):
        f = open(infile)
        lines = f.readlines()
        f.close()
        resultLines = None
        missingLines = None
        signalRegion = ""
        txnames = []
        for il,l in enumerate(lines):
            if not l.strip():
                continue
            if "#Analysis" in l:
                resultLines = True
                continue
            if "Missing topologies" in l:
                missingLines = True
                continue
            if l.startswith("#"):
                continue
            if "----" in l:
                continue
            if "====" in l:
                resultLines = None
                missingLines = None
            if not l.strip():
                continue
            if 'Signal Region' in l or 'Txnames' in l:
                continue
            if 'Missing Topologies' in l or 'Sqrts (TeV)' in l:
                continue
            if resultLines:
                if 'Signal Region' in lines[il+1]:
                    signalRegion = lines[il+1].split(':')[1].strip()
                if 'Txnames' in lines[il+2]:
                    txnames = sorted(lines[il+2].split(':')[1].strip().split(','))
                self.results.append(Output(l.split(),signalRegion,txnames))
            if "Full" in l: break
            if missingLines:
                self.missedTopos.append(MissOutput(l.split()))

    def __eq__(self, other):
        if not len(self.results) == len(other.results):
            return False
        if not len(self.missedTopos) == len(other.missedTopos):
            return False
        for res in self.results:
            if not res in other.results:
                return False
        for miss in self.missedTopos:
            if not miss in other.missedTopos:
                return False
        return True
