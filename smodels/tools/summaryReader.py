"""
.. module:: summaryReader
   :synopsis: Classes to read the summary.txt files.

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""

class Output():
    def __init__(self, l, signalRegion, txnames, allowedDiff=0.):
        self.ana = l[0]
        self.sqrts = l[1]
        self.condViolation = l[2]
        self.tval = l[3]
        self.ul = l[4]
        self.rval = l[5]
        self.txnames = txnames
        self.signalRegion = signalRegion
        self.allowedDiff = allowedDiff

    def __eq__(self, other):
        return fuzzycomp(self,other,self.allowedDiff)

class MissOutput():
    def __init__(self, l,allowedDiff=0.):
        self.sqrts = eval(l[0])
        self.weight = eval(l[1])
        self.element = l[3]
        self.allowedDiff = allowedDiff

    def __eq__(self, other):
        return fuzzycomp(self,other,self.allowedDiff)

        
class AsymmetricOutput():
    def __init__(self, l,allowedDiff=0.):
        self.mother1 = eval(l[0])
        self.mother2 = eval(l[1])
        self.weight = eval(l[2])
        self.allowedDiff = allowedDiff

    def __eq__(self, other):
        return fuzzycomp(self,other,self.allowedDiff)


class Summary():
    """
    Class to access the output given in the summary.txt

    """
    def __init__(self, filename, allowedDiff=0.):
        self.results = []
        self.signalRegions = []
        self.txnames = []
        self.filename = filename
        self.missedTopos = []
        self.outsideTopos = []
        self.longTopos = []
        self.asymmetricTopos = []
        self.allowedDiff = allowedDiff
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
        outsideLines = None
        longLines = None
        asymmetricLines = None
        signalRegion = ""
        txnames = []
        for il,l in enumerate(lines):
            if not l.strip():
                continue
            if "#Analysis" in l:
                resultLines = True
                continue
            elif "Missing topologies" in l:
                missingLines = True
                continue
            elif "Contributions outside" in l:
                outsideLines = True
                continue
            elif "Long cascade decay" in l:
                longLines = True
                continue
            elif "Asymmetric branch decay" in l:
                asymmetricLines = True
                continue
            if l.startswith("#"):
                continue
            if "----" in l:
                continue
            if "====" in l:
                resultLines = None
                missingLines = None
                outsideLines = None
                longLines = None
                asymmetricLines = None                
            if not l.strip():
                continue
            if 'Signal Region' in l or 'Txnames' in l:
                continue
            if (missingLines or outsideLines)  and 'Sqrts (TeV)' in l:
                continue
            if (longLines or asymmetricLines)  and 'Mother1 Mother2' in l:
                continue
            
            if resultLines:
                if 'Signal Region' in lines[il+1]:
                    signalRegion = lines[il+1].split(':')[1].strip()
                if 'Txnames' in lines[il+2]:
                    txnames = sorted(lines[il+2].split(':')[1].strip().split(','))
                self.results.append(Output(l.split(),signalRegion,txnames,self.allowedDiff))

            if missingLines:
                self.missedTopos.append(MissOutput(l.split(),self.allowedDiff))
            if outsideLines:
                self.outsideTopos.append(MissOutput(l.split(),self.allowedDiff))
            if longLines:
                self.longTopos.append(AsymmetricOutput(l.split(),self.allowedDiff))
            if asymmetricLines:
                self.asymmetricTopos.append(AsymmetricOutput(l.split(),self.allowedDiff))

        self.results = sorted(self.results, key=lambda res: [res.ana,res.signalRegion,res.tval])


    def __eq__(self, other):
        if not len(self.results) == len(other.results):            
            return False
        if not len(self.missedTopos) == len(other.missedTopos):
            return False
        if not len(self.outsideTopos) == len(other.outsideTopos):
            return False
        if not len(self.longTopos) == len(other.longTopos):
            return False
        if not len(self.asymmetricTopos) == len(other.asymmetricTopos):
            return False

        if self.results != other.results:
            return False        
        if self.missedTopos != other.missedTopos:
            return False
        if self.outsideTopos != other.outsideTopos:
            return False        
        if self.longTopos != other.longTopos:
            return False        
        if self.asymmetricTopos != other.asymmetricTopos:
            return False
                    
        return True

def fuzzycomp(obj,other,allowedDiff=0.,ignore=[]):
    """
    Compares two objects.
    Accepts numerical differences between floats below allowedDiff.
    """
            
    if isinstance(other, obj.__class__):
        for key,val in obj.__dict__.items():
            if key in ignore: continue
            if not key in other.__dict__:
                return False
            elif val == other.__dict__[key]:
                continue
            else:
                try:
                    val = eval(val)
                    oval = eval(other.__dict__[key])
                    diff = 2.*abs(val-oval)/abs(val+oval)                        
                except (NameError,TypeError):
                    diff = None
                if diff > allowedDiff or diff is None:
                    return False
        return True
    else:
        return False
