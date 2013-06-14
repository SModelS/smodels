#!/usr/bin/env python

"""
.. module:: pyslha2
    :synopsis: missing
    
.. moduleauthor:: missing <email@example.com>
    
"""
    
"""
A simple but flexible parser of SUSY Les Houches Accord (SLHA) model and decay files.

pyslha is a parser/writer module for particle physics SUSY Les Houches Accord
(SLHA) supersymmetric spectrum/decay files, and a collection of scripts which
use the interface, e.g. for conversion to and from the legacy ISAWIG format, or
to plot the mass spectrum and decay chains.

The current release supports SLHA version 1, and as far as we're aware is also
fully compatible with SLHA2: the block structures are read and accessed
completely generically. If you have any problems with SLHA2, please provide an
example input file and we'll investigate.

The plotting script provides output in PDF, EPS and PNG via LaTeX and the TikZ
graphics package, and as LaTeX/TikZ source for direct embedding into documents or
user-tweaking of the generated output.

WARNING: This is a version hacked by Doris Proschofsky. The order of writing out the blocks
has been fixed to satisfy pythia.

TODOs:
 * Split writeSLHA into writeSLHA{Blocks,Decays}
 * Identify HERWIG decay matrix element to use in ISAWIG
 * Handle RPV SUSY in ISAWIG
"""

__author__ = "Andy Buckley <andy.buckley@cern.ch"
__version__ = "1.4.3"


def _autotype(var):
    """Automatically convert strings to numerical types if possible."""
    if type(var) is not str:
        return var
    if var.isdigit() or (var.startswith("-") and var[1:].isdigit()):
        return int(var)
    try:
        f = float(var)
        return f
    except ValueError:
        return var

def _autostr(var, precision=8):
    """Automatically numerical types to the right sort of string."""
    if type(var) is float:
        return ("%." + str(precision) + "e") % var
    return str(var)


class ParseError(Exception):
    "Exception object to be raised when a spectrum file/string is malformed"
    def __init__(self, errmsg):
        self.msg = errmsg
    def __str__(self):
        return self.msg



class Block(object):
    """
    Object representation of any BLOCK elements read from the SLHA file.  Blocks
    have a name, may have an associated Q value, and then a collection of data
    entries, stored as a recursive dictionary. Types in the dictionary are
    numeric (int or float) when a cast from the string in the file has been
    possible.
    """
    def __init__(self, name, q=None):
        self.name = name
        self.entries = {}
        self.q = _autotype(q)

    def add_entry(self, entry):
        #print entry
        nextparent = self.entries
        if len(entry) < 2:
            raise Exception("Block entries must be at least a 2-tuple")
        #print "in", entry
        entry = map(_autotype, entry)
        #print "out", entry
        for e in entry[:-2]:
            if e is not entry[-1]:
                nextparent = nextparent.setdefault(e, {})
        nextparent[entry[-2]] = entry[-1]
        #print self.entries

    def __cmp__(self, other):
        return cmp(self.name, other.name)

    def __str__(self):
        s = self.name
        if self.q is not None:
            s += " (Q=%s)" % self.q
        s += "\n"
        s += str(self.entries)
        return s

    def __repr__(self):
        return self.__str__()


class Decay(object):
    """
    Object representing a decay entry on a particle decribed by the SLHA file.
    'Decay' objects are not a direct representation of a DECAY block in an SLHA
    file... that role, somewhat confusingly, is taken by the Particle class.

    Decay objects have three properties: a branching ratio, br, an nda number
    (number of daughters == len(ids)), and a tuple of PDG PIDs to which the
    decay occurs. The PDG ID of the particle whose decay this represents may
    also be stored, but this is normally known via the Particle in which the
    decay is stored.
    """
    def __init__(self, br, nda, ids, parentid=None):
        self.parentid = parentid
        self.br = br
        self.nda = nda
        self.ids = ids
        assert(self.nda == len(self.ids))

    def __cmp__(self, other):
        return cmp(other.br, self.br)

    def __str__(self):
        return "%.8e %s" % (self.br, self.ids)

    def __repr__(self):
        return self.__str__()


class Particle(object):
    """
    Representation of a single, specific particle, decay block from an SLHA
    file.  These objects are not themselves called 'Decay', since that concept
    applies more naturally to the various decays found inside this
    object. Particle classes store the PDG ID (pid) of the particle being
    represented, and optionally the mass (mass) and total decay width
    (totalwidth) of that particle in the SLHA scenario. Masses may also be found
    via the MASS block, from which the Particle.mass property is filled, if at
    all. They also store a list of Decay objects (decays) which are probably the
    item of most interest.
    """
    def __init__(self, pid, totalwidth=None, mass=None):
        self.pid = pid
        self.totalwidth = totalwidth
        self.mass = mass
        self.decays = []

    def add_decay(self, br, nda, ids):
        self.decays.append(Decay(br, nda, ids))
        self.decays.sort()

    def __cmp__(self, other):
        if abs(self.pid) == abs(other.pid):
            return cmp(self.pid, other.pid)
        return cmp(abs(self.pid), abs(other.pid))

    def __str__(self):
        s = str(self.pid)
        if self.mass is not None:
            s += " : mass = %.8e GeV" % self.mass
        if self.totalwidth is not None:
            s += " : total width = %.8e GeV" % self.totalwidth
        for d in self.decays:
            if d.br > 0.0:
                s += "\n  %s" % d
        return s

    def __repr__(self):
        return self.__str__()





def readSLHA(spcstr, ignorenobr=False):
    """
    Read an SLHA definition from a string, returning dictionaries of blocks and
    decays.

    If the ignorenobr parameter is True, do not store decay entries with a
    branching ratio of zero.
    """
    blocks = {}
    decays = {}
    #
    listBlocks=[]
    listDecays=[]
    #
    import re
    currentblock = None
    currentdecay = None
    for line in spcstr.splitlines():
        ## Handle (ignore) comment lines
        if line.startswith("#"):
            continue
        if "#" in line:
            line = line[:line.index("#")]

        ## Handle BLOCK/DECAY start lines
        if line.upper().startswith("BLOCK"):
            #print line
            match = re.match(r"BLOCK\s+(\w+)(\s+Q\s*=\s*.+)?", line.upper())
            if not match:
                continue
            blockname = match.group(1)
            listBlocks.append(blockname)
            qstr = match.group(2)
            if qstr is not None:
                qstr = qstr[qstr.find("=")+1:].strip()
            currentblock = blockname
            currentdecay = None
            blocks[blockname] = Block(blockname, q=qstr)
        elif line.upper().startswith("DECAY"):
            match = re.match(r"DECAY\s+(\d+)\s+([\d\.E+-]+).*", line.upper())
            if not match:
                continue
            pdgid = int(match.group(1))
            listDecays.append(pdgid)
            width = float(match.group(2))
            currentblock = "DECAY"
            currentdecay = pdgid
            decays[pdgid] = Particle(pdgid, width)
        else:
            ## In-block line
            if currentblock is not None:
                items = line.split()
                if len(items) < 1:
                    continue
                if currentblock != "DECAY":
                    if len(items) < 2:
                        ## Treat the ALPHA block differently
                        blocks[currentblock].value = _autotype(items[0])
                        blocks[currentblock].entries = _autotype(items[0])
                    else:
                        blocks[currentblock].add_entry(items)
                else:
                    br = float(items[0])
                    nda = int(items[1])
                    ids = map(int, items[2:])
                    if br > 0.0 or not ignorenobr:
                        decays[currentdecay].add_decay(br, nda, ids)

    ## Try to populate Particle masses from the MASS block
    # print blocks.keys()
    try:
        for pid in blocks["MASS"].entries.keys():
            if decays.has_key(pid):
                decays[pid].mass = blocks["MASS"].entries[pid]
    except:
        raise ParseError("No MASS block found: cannot populate particle masses")

#    print listBlocks, listDecays
    return blocks, decays, listBlocks, listDecays




# TODO: Split writeSLHA into writeSLHA{Blocks,Decays}


#def writeSLHA(blocks, decays, listBlocks, listDecays, ignorenobr=False, precision=8):
def writeSLHA(f, ignorenobr=False, precision=8):
    """
    Return an SLHA definition as a string, from the supplied blocks and decays dicts.
    """
    fmte = "%." + str(precision) + "e"

    sep = "   "
    out = ""
    blocks=f[0]
    decays=f[1]
    listBlocks=f[2]
    listDecays=f[3]
#    listBlocks=["DCINFO", "SPINFO", "MODSEL", "SMINPUTS", "MINPAR", "EXTPAR", "MASS", "NMIX", "UMIX", "VMIX", "STOPMIX", "SBOTMIX", "STAUMIX", "ALPHA", "HMIX", "GAUGE", "AU", "AD", "AE", "Yu", "Yd", "Ye", "MSOFT" ]
#    listDecays=["", "", "", "", "", "", "", "", ""]
    def dict_hier_strs(d, s=""):
        if type(d) is dict: 
            for k, v in sorted(d.iteritems()):	### iteritems(): (key, value)-list of dict
                for s2 in dict_hier_strs(v, s + sep + _autostr(k)):
                    yield s2
        else:
            yield s + sep + _autostr(d)
    ## Blocks
    for bname in listBlocks:
        namestr = bname
        if blocks[bname].q is not None:
            namestr += (" Q=   " + fmte) % float(blocks[bname].q)
        out += "BLOCK %s\n" % namestr
	if bname == "SPINFO":
	    for k, v in blocks[bname].entries.iteritems():
		count = 6 - len(str(k))
		space = ""
		for x in range(count):
		     space += " "
		out += space + str(k) + sep + v + "\n"
	    out += "\n"
	else:
	    for s in dict_hier_strs(blocks[bname].entries):
        	out += sep + s + "\n"
       	    out += "\n"
        
#    for bname, b in sorted(blocks.iteritems()):
#        namestr = b.name
#        if b.q is not None:
#            namestr += (" Q= " + fmte) % float(b.q)
#        out += "BLOCK %s\n" % namestr
#        for s in dict_hier_strs(b.entries):
#            out += sep + s + "\n"
#        out += "\n"
    ## Decays
    for dname in listDecays:
        pid = dname
        out += ("DECAY %d " + fmte + "\n") % (decays[pid].pid, decays[pid].totalwidth or -1)
        for d in sorted(decays[pid].decays):
            if d.br > 0.0 or not ignorenobr:
                products_str = "   ".join(map(str, d.ids))
                out += sep + fmte % d.br + sep + "%d" % len(d.ids) + sep + products_str + "\n"
        out += "\n"
    return out

#    for pid, particle in sorted(decays.iteritems()):
#        out += ("DECAY %d " + fmte + "\n") % (particle.pid, particle.totalwidth or -1)
#        for d in sorted(particle.decays):
#           if d.br > 0.0 or not ignorenobr:
#                products_str = "   ".join(map(str, d.ids))
#                out += sep + fmte % d.br + sep + "%d" % len(d.ids) + sep + products_str + "\n"
#        out += "\n"
#    return out



###############################################################################
## PDG <-> HERWIG particle ID code translations for ISAWIG handling

## Static array of HERWIG IDHW codes mapped to PDG MC ID codes, based on
## http://www.hep.phy.cam.ac.uk/~richardn/HERWIG/ISAWIG/susycodes.html
## + the IDPDG array and section 4.13 of the HERWIG manual.
_HERWIGID2PDGID = {}
_HERWIGID2PDGID[7]   = -1
_HERWIGID2PDGID[8]   = -2
_HERWIGID2PDGID[9]   = -3
_HERWIGID2PDGID[10]  = -4
_HERWIGID2PDGID[11]  = -5
_HERWIGID2PDGID[12]  = -6
_HERWIGID2PDGID[13]  =  21
_HERWIGID2PDGID[59]  =  22
_HERWIGID2PDGID[121] =  11
_HERWIGID2PDGID[122] =  12
_HERWIGID2PDGID[123] =  13
_HERWIGID2PDGID[124] =  14
_HERWIGID2PDGID[125] =  15
_HERWIGID2PDGID[126] =  16
_HERWIGID2PDGID[127] = -11
_HERWIGID2PDGID[128] = -12
_HERWIGID2PDGID[129] = -13
_HERWIGID2PDGID[130] = -14
_HERWIGID2PDGID[131] = -15
_HERWIGID2PDGID[132] = -16
_HERWIGID2PDGID[198] =  24 # W+
_HERWIGID2PDGID[199] = -24 # W-
_HERWIGID2PDGID[200] =  23 # Z0
_HERWIGID2PDGID[201] =  25 ## SM HIGGS
_HERWIGID2PDGID[203] =  25 ## HIGGSL0 (== PDG standard in this direction)
_HERWIGID2PDGID[204] =  35 ## HIGGSH0
_HERWIGID2PDGID[205] =  36 ## HIGGSA0
_HERWIGID2PDGID[206] =  37 ## HIGGS+
_HERWIGID2PDGID[207] = -37 ## HIGGS-
_HERWIGID2PDGID[401] =  1000001 ## SSDLBR
_HERWIGID2PDGID[407] = -1000001 ## SSDLBR
_HERWIGID2PDGID[402] =  1000002 ## SSULBR
_HERWIGID2PDGID[408] = -1000002 ## SSUL
_HERWIGID2PDGID[403] =  1000003 ## SSSLBR
_HERWIGID2PDGID[409] = -1000003 ## SSSL
_HERWIGID2PDGID[404] =  1000004 ## SSCLBR
_HERWIGID2PDGID[410] = -1000004 ## SSCL
_HERWIGID2PDGID[405] =  1000005 ## SSB1BR
_HERWIGID2PDGID[411] = -1000005 ## SSB1
_HERWIGID2PDGID[406] =  1000006 ## SST1BR
_HERWIGID2PDGID[412] = -1000006 ## SST1
_HERWIGID2PDGID[413] =  2000001 ## SSDR
_HERWIGID2PDGID[419] = -2000001 ## SSDRBR
_HERWIGID2PDGID[414] =  2000002 ## SSUR
_HERWIGID2PDGID[420] = -2000002 ## SSURBR
_HERWIGID2PDGID[415] =  2000003 ## SSSR
_HERWIGID2PDGID[421] = -2000003 ## SSSRBR
_HERWIGID2PDGID[416] =  2000004 ## SSCR
_HERWIGID2PDGID[422] = -2000004 ## SSCRBR
_HERWIGID2PDGID[417] =  2000005 ## SSB2
_HERWIGID2PDGID[423] = -2000005 ## SSB2BR
_HERWIGID2PDGID[418] =  2000006 ## SST2
_HERWIGID2PDGID[424] = -2000006 ## SST2BR
_HERWIGID2PDGID[425] =  1000011 ## SSEL-
_HERWIGID2PDGID[431] = -1000011 ## SSEL+
_HERWIGID2PDGID[426] =  1000012 ## SSNUEL
_HERWIGID2PDGID[432] = -1000012 ## SSNUELBR
_HERWIGID2PDGID[427] =  1000013 ## SSMUL-
_HERWIGID2PDGID[433] = -1000013 ## SSMUL+
_HERWIGID2PDGID[428] =  1000014 ## SSNUMUL
_HERWIGID2PDGID[434] = -1000014 ## SSNUMLBR
_HERWIGID2PDGID[429] =  1000015 ## SSTAU1-
_HERWIGID2PDGID[435] = -1000015 ## SSTAU1+
_HERWIGID2PDGID[430] =  1000016 ## SSNUTL
_HERWIGID2PDGID[436] = -1000016 ## SSNUTLBR
_HERWIGID2PDGID[437] =  2000011 ## SSEL-
_HERWIGID2PDGID[443] = -2000011 ## SSEL+
_HERWIGID2PDGID[438] =  2000012 ## SSNUEL
_HERWIGID2PDGID[444] = -2000012 ## SSNUELBR
_HERWIGID2PDGID[439] =  2000013 ## SSMUL-
_HERWIGID2PDGID[445] = -2000013 ## SSMUL+
_HERWIGID2PDGID[440] =  2000014 ## SSNUMUL
_HERWIGID2PDGID[446] = -2000014 ## SSNUMLBR
_HERWIGID2PDGID[441] =  2000015 ## SSTAU1-
_HERWIGID2PDGID[447] = -2000015 ## SSTAU1+
_HERWIGID2PDGID[442] =  2000016 ## SSNUTL
_HERWIGID2PDGID[448] = -2000016 ## SSNUTLBR
_HERWIGID2PDGID[449] =  1000021 ## GLUINO
_HERWIGID2PDGID[450] =  1000022 ## NTLINO1
_HERWIGID2PDGID[451] =  1000023 ## NTLINO2
_HERWIGID2PDGID[452] =  1000025 ## NTLINO3
_HERWIGID2PDGID[453] =  1000035 ## NTLINO4
_HERWIGID2PDGID[454] =  1000024 ## CHGINO1+
_HERWIGID2PDGID[456] = -1000024 ## CHGINO1-
_HERWIGID2PDGID[455] =  1000037 ## CHGINO2+
_HERWIGID2PDGID[457] = -1000037 ## CHGINO2-
_HERWIGID2PDGID[458] =  1000039 ## GRAVTINO

def herwigid2pdgid(hwid):
    """
    Convert a particle ID code in the HERWIG internal IDHW format (as used by
    ISAWIG) into its equivalent in the standard PDG ID code definition.
    """
    return _HERWIGID2PDGID.get(hwid, hwid)


## PDG MC ID codes mapped to HERWIG IDHW codes, based on
## http://www.hep.phy.cam.ac.uk/~richardn/HERWIG/ISAWIG/susycodes.html
## + the IDPDG array and section 4.13 of the HERWIG manual.
_PDGID2HERWIGID = {}
_PDGID2HERWIGID[      -1] = 7
_PDGID2HERWIGID[      -2] = 8
_PDGID2HERWIGID[      -3] = 9
_PDGID2HERWIGID[      -4] = 10
_PDGID2HERWIGID[      -5] = 11
_PDGID2HERWIGID[      -6] = 12
_PDGID2HERWIGID[      21] = 13
_PDGID2HERWIGID[      22] = 59
_PDGID2HERWIGID[      11] = 121
_PDGID2HERWIGID[      12] = 122
_PDGID2HERWIGID[      13] = 123
_PDGID2HERWIGID[      14] = 124
_PDGID2HERWIGID[      15] = 125
_PDGID2HERWIGID[      16] = 126
_PDGID2HERWIGID[     -11] = 127
_PDGID2HERWIGID[     -12] = 128
_PDGID2HERWIGID[     -13] = 129
_PDGID2HERWIGID[     -14] = 130
_PDGID2HERWIGID[     -15] = 131
_PDGID2HERWIGID[     -16] = 132
_PDGID2HERWIGID[      24] = 198 ## W+
_PDGID2HERWIGID[     -24] = 199 ## W-
_PDGID2HERWIGID[      23] = 200 ## Z
_PDGID2HERWIGID[      25] = 203 ## HIGGSL0 (added for PDG standard -> HERWIG IDHW) # TODO: should be 201?
_PDGID2HERWIGID[      26] = 203 ## HIGGSL0
_PDGID2HERWIGID[      35] = 204 ## HIGGSH0
_PDGID2HERWIGID[      36] = 205 ## HIGGSA0
_PDGID2HERWIGID[      37] = 206 ## HIGGS+
_PDGID2HERWIGID[     -37] = 207 ## HIGGS-
_PDGID2HERWIGID[ 1000001] = 401 ## SSDLBR
_PDGID2HERWIGID[-1000001] = 407 ## SSDLBR
_PDGID2HERWIGID[ 1000002] = 402 ## SSULBR
_PDGID2HERWIGID[-1000002] = 408 ## SSUL
_PDGID2HERWIGID[ 1000003] = 403 ## SSSLBR
_PDGID2HERWIGID[-1000003] = 409 ## SSSL
_PDGID2HERWIGID[ 1000004] = 404 ## SSCLBR
_PDGID2HERWIGID[-1000004] = 410 ## SSCL
_PDGID2HERWIGID[ 1000005] = 405 ## SSB1BR
_PDGID2HERWIGID[-1000005] = 411 ## SSB1
_PDGID2HERWIGID[ 1000006] = 406 ## SST1BR
_PDGID2HERWIGID[-1000006] = 412 ## SST1
_PDGID2HERWIGID[ 2000001] = 413 ## SSDR
_PDGID2HERWIGID[-2000001] = 419 ## SSDRBR
_PDGID2HERWIGID[ 2000002] = 414 ## SSUR
_PDGID2HERWIGID[-2000002] = 420 ## SSURBR
_PDGID2HERWIGID[ 2000003] = 415 ## SSSR
_PDGID2HERWIGID[-2000003] = 421 ## SSSRBR
_PDGID2HERWIGID[ 2000004] = 416 ## SSCR
_PDGID2HERWIGID[-2000004] = 422 ## SSCRBR
_PDGID2HERWIGID[ 2000005] = 417 ## SSB2
_PDGID2HERWIGID[-2000005] = 423 ## SSB2BR
_PDGID2HERWIGID[ 2000006] = 418 ## SST2
_PDGID2HERWIGID[-2000006] = 424 ## SST2BR
_PDGID2HERWIGID[ 1000011] = 425 ## SSEL-
_PDGID2HERWIGID[-1000011] = 431 ## SSEL+
_PDGID2HERWIGID[ 1000012] = 426 ## SSNUEL
_PDGID2HERWIGID[-1000012] = 432 ## SSNUELBR
_PDGID2HERWIGID[ 1000013] = 427 ## SSMUL-
_PDGID2HERWIGID[-1000013] = 433 ## SSMUL+
_PDGID2HERWIGID[ 1000014] = 428 ## SSNUMUL
_PDGID2HERWIGID[-1000014] = 434 ## SSNUMLBR
_PDGID2HERWIGID[ 1000015] = 429 ## SSTAU1-
_PDGID2HERWIGID[-1000015] = 435 ## SSTAU1+
_PDGID2HERWIGID[ 1000016] = 430 ## SSNUTL
_PDGID2HERWIGID[-1000016] = 436 ## SSNUTLBR
_PDGID2HERWIGID[ 2000011] = 437 ## SSEL-
_PDGID2HERWIGID[-2000011] = 443 ## SSEL+
_PDGID2HERWIGID[ 2000012] = 438 ## SSNUEL
_PDGID2HERWIGID[-2000012] = 444 ## SSNUELBR
_PDGID2HERWIGID[ 2000013] = 439 ## SSMUL-
_PDGID2HERWIGID[-2000013] = 445 ## SSMUL+
_PDGID2HERWIGID[ 2000014] = 440 ## SSNUMUL
_PDGID2HERWIGID[-2000014] = 446 ## SSNUMLBR
_PDGID2HERWIGID[ 2000015] = 441 ## SSTAU1-
_PDGID2HERWIGID[-2000015] = 447 ## SSTAU1+
_PDGID2HERWIGID[ 2000016] = 442 ## SSNUTL
_PDGID2HERWIGID[-2000016] = 448 ## SSNUTLBR
_PDGID2HERWIGID[ 1000021] = 449 ## GLUINO
_PDGID2HERWIGID[ 1000022] = 450 ## NTLINO1
_PDGID2HERWIGID[ 1000023] = 451 ## NTLINO2
_PDGID2HERWIGID[ 1000025] = 452 ## NTLINO3
_PDGID2HERWIGID[ 1000035] = 453 ## NTLINO4
_PDGID2HERWIGID[ 1000024] = 454 ## CHGINO1+
_PDGID2HERWIGID[-1000024] = 456 ## CHGINO1-
_PDGID2HERWIGID[ 1000037] = 455 ## CHGINO2+
_PDGID2HERWIGID[-1000037] = 457 ## CHGINO2-
_PDGID2HERWIGID[ 1000039] = 458 ## GRAVTINO

def pdgid2herwigid(pdgid):
    """
    Convert a particle ID code in the standard PDG ID code definition into
    its equivalent in the HERWIG internal IDHW format (as used by ISAWIG).
    """
    return _PDGID2HERWIGID.get(pdgid, pdgid)


###############################################################################
## ISAWIG format reading/writing


def readISAWIG(isastr, ignorenobr=False):
    """
    Read a spectrum definition from a string in the ISAWIG format, returning
    dictionaries of blocks and decays. While this is not an SLHA format, it is
    informally supported as a useful mechanism for converting ISAWIG spectra to
    SLHA.

    ISAWIG parsing based on the HERWIG SUSY specification format, from
    http://www.hep.phy.cam.ac.uk/~richardn/HERWIG/ISAWIG/file.html

    If the ignorenobr parameter is True, do not store decay entries with a
    branching ratio of zero.
    """

    blocks = {}
    decays = {}
    LINES = isastr.splitlines()

    def getnextvalidline():
        while LINES:
            s = LINES.pop(0).strip()
            ## Return None if EOF reached
            if len(s) == 0:
                continue
            ## Strip comments
            if "#" in s:
                s = s[:s.index("#")].strip()
            ## Return if non-empty
            if len(s) > 0:
                return s

    def getnextvalidlineitems():
        return map(_autotype, getnextvalidline().split())

    ## Populate MASS block and create decaying particle objects
    masses = Block("MASS")
    numentries = int(getnextvalidline())
    for i in xrange(numentries):
        hwid, mass, lifetime = getnextvalidlineitems()
        width = 1.0/(lifetime * 1.51926778e24) ## width in GeV == hbar/lifetime in seconds
        pdgid = herwigid2pdgid(hwid)
        masses.add_entry((pdgid, mass))
        decays[pdgid] = Particle(pdgid, width, mass)
        #print pdgid, mass, width
    blocks["MASS"] = masses

    ## Populate decays
    for n in xrange(numentries):
        numdecays = int(getnextvalidline())
        for d in xrange(numdecays):
            #print n, numentries-1, d, numdecays-1
            decayitems = getnextvalidlineitems()
            hwid = decayitems[0]
            pdgid = herwigid2pdgid(hwid)
            br = decayitems[1]
            nme = decayitems[2]
            daughter_hwids = decayitems[3:]
            daughter_pdgids = []
            for hw in daughter_hwids:
                if hw != 0:
                    daughter_pdgids.append(herwigid2pdgid(hw))
            if not decays.has_key(pdgid):
                #print "Decay for unlisted particle %d, %d" % (hwid, pdgid)
                decays[pdgid] = Particle(pdgid)
            decays[pdgid].add_decay(br, len(daughter_pdgids), daughter_pdgids)


    ## Now the SUSY parameters
    TANB, ALPHAH = getnextvalidlineitems()
    blocks["MINPAR"] = Block("MINPAR")
    blocks["MINPAR"].add_entry((3, TANB))
    blocks["ALPHA"] = Block("ALPHA")
    blocks["ALPHA"].entries = ALPHAH
    #
    ## Neutralino mixing matrix
    blocks["NMIX"] = Block("NMIX")
    for i in xrange(1, 5):
        nmix_i = getnextvalidlineitems()
        for j, v in enumerate(nmix_i):
            blocks["NMIX"].add_entry((i, j+1, v))
    #
    ## Chargino mixing matrices V and U
    blocks["VMIX"] = Block("VMIX")
    vmix = getnextvalidlineitems()
    blocks["VMIX"].add_entry((1, 1, vmix[0]))
    blocks["VMIX"].add_entry((1, 2, vmix[1]))
    blocks["VMIX"].add_entry((2, 1, vmix[2]))
    blocks["VMIX"].add_entry((2, 2, vmix[3]))
    blocks["UMIX"] = Block("UMIX")
    umix = getnextvalidlineitems()
    blocks["UMIX"].add_entry((1, 1, umix[0]))
    blocks["UMIX"].add_entry((1, 2, umix[1]))
    blocks["UMIX"].add_entry((2, 1, umix[2]))
    blocks["UMIX"].add_entry((2, 2, umix[3]))
    #
    THETAT, THETAB, THETAL = getnextvalidlineitems()
    import math
    blocks["STOPMIX"] = Block("STOPMIX")
    blocks["STOPMIX"].add_entry((1, 1,  math.cos(THETAT)))
    blocks["STOPMIX"].add_entry((1, 2, -math.sin(THETAT)))
    blocks["STOPMIX"].add_entry((2, 1,  math.sin(THETAT)))
    blocks["STOPMIX"].add_entry((2, 2,  math.cos(THETAT)))
    blocks["SBOTMIX"] = Block("SBOTMIX")
    blocks["SBOTMIX"].add_entry((1, 1,  math.cos(THETAB)))
    blocks["SBOTMIX"].add_entry((1, 2, -math.sin(THETAB)))
    blocks["SBOTMIX"].add_entry((2, 1,  math.sin(THETAB)))
    blocks["SBOTMIX"].add_entry((2, 2,  math.cos(THETAB)))
    blocks["STAUMIX"] = Block("STAUMIX")
    blocks["STAUMIX"].add_entry((1, 1,  math.cos(THETAL)))
    blocks["STAUMIX"].add_entry((1, 2, -math.sin(THETAL)))
    blocks["STAUMIX"].add_entry((2, 1,  math.sin(THETAL)))
    blocks["STAUMIX"].add_entry((2, 2,  math.cos(THETAL)))
    #
    ATSS, ABSS, ALSS = getnextvalidlineitems()
    blocks["AU"] = Block("AU")
    blocks["AU"].add_entry((3, 3, ATSS))
    blocks["AD"] = Block("AD")
    blocks["AD"].add_entry((3, 3, ABSS))
    blocks["AE"] = Block("AE")
    blocks["AE"].add_entry((3, 3, ALSS))
    #
    MUSS = getnextvalidlineitems()[0]
    blocks["MINPAR"].add_entry((4, MUSS))
    #

    # TODO: Parse RPV boolean and couplings into SLHA2 blocks

    return blocks, decays


def writeISAJET(blocks, decays, outname, ignorenobr=False, precision=8):
    """
    Return a SUSY spectrum definition in the format required for input by ISAJET,
    as a string, from the supplied blocks and decays dicts.

    The outname parameter specifies the desired output filename from ISAJET: this
    will appear in the first line of the return value.

    If the ignorenobr parameter is True, do not write decay entries with a
    branching ratio of zero.
    """
    fmte = "%." + str(precision) + "e"

    masses = blocks["MASS"].entries

    ## Init output string
    out = ""

    ## First line is the output name
    out += "'%s'" % outname + "\n"

    ## Next the top mass
    out += fmte % masses[6] + "\n"

    ## Next the top mass
    out += fmte % masses[6] + "\n"

    ## mSUGRA parameters (one line)
    # e.g. 1273.78,713.286,804.721,4.82337

    ## Masses and trilinear couplings (3 lines)
    # e.g. 1163.14,1114.15,1118.99,374.664,209.593
    # e.g. 1069.54,1112.7,919.908,374.556,209.381,-972.817,-326.745,-406.494
    # e.g. 1163.14,1114.15,1118.99,374.712,210.328

    ## RPV couplings (?? lines)
    # e.g. 232.615,445.477

    ## Etc ???!!!
    # e.g. /
    # e.g. n
    # e.g. y
    # e.g. y
    # e.g. 0.047441 3.80202e-23 0 0 0 2.17356e-22 0 0 5.23773e-09
    # e.g. y
    # e.g. 3.35297e-25 0 0 0 7.34125e-24 0 0 0 3.17951e-22 8.07984e-12 0 0 0 1.76906e-10 0 0 0 7.66184e-09 0 0 0 0 0 0 0 0 0
    # e.g. n
    # e.g. 'susy_RPV_stau_BC1scan_m560_tanb05.txt'

    return out


def writeISAWIG(blocks, decays, ignorenobr=False, precision=8):
    """
    Return a SUSY spectrum definition in the format produced by ISAWIG for inut to HERWIG
    as a string, from the supplied SLHA blocks and decays dicts.

    ISAWIG parsing based on the HERWIG SUSY specification format, from
    http://www.hep.phy.cam.ac.uk/~richardn/HERWIG/ISAWIG/file.html

    If the ignorenobr parameter is True, do not write decay entries with a
    branching ratio of zero.
    """
    fmte = "%." + str(precision) + "e"

    masses = blocks["MASS"].entries

    ## Init output string
    out = ""

    ## First write out masses section:
    ##   Number of SUSY + top particles
    ##   IDHW, RMASS(IDHW), RLTIM(IDHW)
    ##   repeated for each particle
    ## IDHW is the HERWIG identity code.
    ## RMASS and RTLIM are the mass in GeV, and lifetime in seconds respectively.
    massout = ""
    for pid in masses.keys():
        lifetime = -1
        try:
            width = decays[pid].totalwidth
            if width and width > 0:
                lifetime = 1.0/(width * 1.51926778e24) ## lifetime in seconds == hbar/width in GeV
        except:
            pass
        massout += ("%d " + fmte + " " + fmte + "\n") % (pdgid2herwigid(pid), masses[pid], lifetime)
    out += "%d\n" % massout.count("\n")
    out += massout

    assert(len(masses) == len(decays))

    ## Next each particles decay modes together with their branching ratios and matrix element codes
    ##   Number of decay modes for a given particle (IDK)
    ##     IDK(*), BRFRAC(*), NME(*) & IDKPRD(1-5,*)
    ##     repeated for each mode.
    ##   Repeated for each particle.
    ## IDK is the HERWIG code for the decaying particle, BRFRAC is the branching ratio of
    ## the decay mode. NME is a code for the matrix element to be used, either from the
    ## SUSY elements or the main HERWIG MEs. IDKPRD are the HERWIG identity codes of the decay products.
    for i, pid in enumerate(decays.keys()):
        # if not decays.has_key(pid):
        #     continue
        hwid = pdgid2herwigid(pid)
        decayout = ""
        #decayout += "@@@@ %d %d %d\n" % (i, pid, hwid)
        for i_d, d in enumerate(decays[pid].decays):
            ## Skip decay if it has no branching ratio
            if ignorenobr and d.br == 0:
                continue

            ## Identify decay matrix element to use
            ## From std HW docs, or from this pair:
            ## Two new matrix element codes have been added for these new decays:
            ##    NME =	200 	3 body top quark via charged Higgs
            ##    	300 	3 body R-parity violating gaugino and gluino decays
            nme = 0
            # TODO: Get correct condition for using ME 100... this guessed from some ISAWIG output
            if abs(pid) in (6, 12):
                nme = 100
            ## Extra SUSY MEs
            if len(d.ids) == 3:
                # TODO: How to determine the conditions for using 200 and 300 MEs? Enumeration of affected decays?
                pass
            decayout += "%d " + fmte + " %d " % (hwid, d.br, nme)

            def is_quark(pid):
                return (abs(pid) in range(1, 7))

            def is_lepton(pid):
                return (abs(pid) in range(11, 17))

            def is_squark(pid):
                if abs(pid) in range(1000001, 1000007):
                    return True
                if abs(pid) in range(2000001, 2000007):
                    return True
                return False

            def is_slepton(pid):
                if abs(pid) in range(1000011, 1000017):
                    return True
                if abs(pid) in range(2000011, 2000016, 2):
                    return True
                return False

            def is_gaugino(pid):
                if abs(pid) in range(1000022, 1000026):
                    return True
                if abs(pid) in (1000035, 1000037):
                    return True
                return False

            def is_susy(pid):
                return (is_squark(pid) or is_slepton(pid) or is_gaugino(pid) or pid == 1000021)

            absids = map(abs, d.ids)

            ## Order decay products as required by HERWIG
            ## Top
            if abs(pid) == 6:
                def cmp_bottomlast(a, b):
                    """Comparison function which always puts b/bbar last"""
                    if abs(a) == 5:
                        return True
                    if abs(b) == 5:
                        return False
                    return cmp(a, b)
                if len(absids) == 2:
                    ## 2 body mode, to Higgs: Higgs; Bottom
                    if (25 in absids or 26 in absids) and 5 in absids:
                        d.ids = sorted(d.ids, key=cmp_bottomlast)
                elif len(absids) == 3:
                    ## 3 body mode, via charged Higgs/W: quarks or leptons from W/Higgs; Bottom
                    if 37 in absids or 23 in absids:
                        d.ids = sorted(d.ids, key=cmp_bottomlast)
            ## Gluino
            elif abs(pid) == 1000021:
                if len(absids) == 2:
                    ## 2 body mode
                    ## without gluon: any order
                    ## with gluon: gluon; colour neutral
                    if 21 in absids:
                        def cmp_gluonfirst(a, b):
                            """Comparison function which always puts gluon first"""
                            if a == 21:
                                return False
                            if b == 21:
                                return True
                            return cmp(a, b)
                        d.ids = sorted(d.ids, key=cmp_gluonfirst)
                elif len(absids) == 3:
                    ## 3-body modes, R-parity conserved: colour neutral; q or qbar
                    def cmp_quarkslast(a, b):
                        """Comparison function which always puts quarks last"""
                        if is_quark(a):
                            return True
                        if is_quark(b):
                            return False
                        return cmp(a, b)
                    d.ids = sorted(d.ids, key=cmp_quarkslast)
            ## Squark/Slepton
            elif is_squark(pid) or is_slepton(pid):
                def cmp_susy_quark_lepton(a, b):
                    if is_susy(a):
                        return False
                    if is_susy(b):
                        return True
                    if is_quark(a):
                        return False
                    if is_quark(b):
                        return True
                    return cmp(a, b)
                ##   2 body modes: Gaugino/Gluino with Quark/Lepton     Gaugino      quark
                ##                                                      Gluino       lepton
                ##   3 body modes: Weak                                 sparticle    particles from W decay
                ## Squark
                ##   2 body modes: Lepton Number Violated               quark     lepton
                ##                 Baryon number violated               quark     quark
                ## Slepton
                ##   2 body modes: Lepton Number Violated               q or qbar
                d.ids = sorted(d.ids, key=cmp_bottomlast)
            ## Higgs
            elif pid in (25, 26):
                # TODO: Includes SUSY Higgses?
                ## Higgs
                ##   2 body modes: (s)quark-(s)qbar                     (s)q or (s)qbar
                ##                 (s)lepton-(s)lepton                  (s)l or (s)lbar
                ##   3 body modes:                                      colour neutral       q or qbar
                if len(absids) == 3:
                    def cmp_quarkslast(a, b):
                        """Comparison function which always puts quarks last"""
                        if is_quark(a):
                            return True
                        if is_quark(b):
                            return False
                        return cmp(a, b)
                    d.ids = sorted(d.ids, key=cmp_quarkslast)
            elif is_gaugino(pid):
                # TODO: Is there actually anything to do here?
                ## Gaugino
                ##   2 body modes: Squark-quark                         q or sq
                ##                 Slepton-lepton                       l or sl
                ##
                ##   3 body modes: R-parity conserved                   colour neutral       q or qbar
                ##                                                                           l or lbar
                if len(absids) == 3:
                    def cmp_quarkslast(a, b):
                        """Comparison function which always puts quarks last"""
                        if is_quark(a):
                            return True
                        if is_quark(b):
                            return False
                        return cmp(a, b)
                    d.ids = sorted(d.ids, key=cmp_quarkslast)

            # TODO: Gaugino/Gluino
            ##   3 body modes:  R-parity violating:   Particles in the same order as the R-parity violating superpotential

            ## Pad out IDs list with zeros
            ids = [0,0,0,0,0]
            for i, pid in enumerate(d.ids):
                ids[i] = pid
            ids = map(str, ids)
            decayout += " ".join(ids) + "\n"
        decayout = "%d\n" % decayout.count("\n") + decayout
        out += decayout

    ## Now the SUSY parameters
    ## TANB, ALPHAH:
    out += (fmte + " " + fmte + "\n") % (blocks["MINPAR"].entries[3], blocks["ALPHA"].entries)
    ## Neutralino mixing matrix
    nmix = blocks["NMIX"].entries
    for i in xrange(1, 5):
        out += (fmte + " " + fmte + " " + fmte + " " + fmte + "\n") % (nmix[i][1], nmix[i][2], nmix[i][3], nmix[i][4])
    ## Chargino mixing matrices V and U
    vmix = blocks["VMIX"].entries
    out += (fmte + " " + fmte + " " + fmte + " " + fmte + "\n") % (vmix[1][1], vmix[1][2], vmix[2][1], vmix[2][2])
    umix = blocks["UMIX"].entries
    out += (fmte + " " + fmte + " " + fmte + " " + fmte + "\n") % (umix[1][1], umix[1][2], umix[2][1], umix[2][2])
    ## THETAT,THETAB,THETAL
    import math
    out += (fmte + " " + fmte + " " + fmte + " " + "\n") % (math.acos(blocks["STOPMIX"].entries[1][1]),
                                                            math.acos(blocks["SBOTMIX"].entries[1][1]),
                                                            math.acos(blocks["STAUMIX"].entries[1][1]))
    ## ATSS,ABSS,ALSS
    out += (fmte + " " + fmte + " " + fmte + " " + "\n") % (blocks["AU"].entries[3][3],
                                                            blocks["AD"].entries[3][3],
                                                            blocks["AE"].entries[3][3])
    ## MUSS == sign(mu)
    out += "%f\n" % blocks["MINPAR"].entries[4]

    ## RPV SUSY
    isRPV = False
    out += "%d\n" % isRPV
    # TODO: Write RPV couplings if RPV is True (lambda1,2,3; 27 params in each, sci format.
    # TODO: Get the index orderings right
    # if isRPV: ...

    return out


###############################################################################
## File-level functions


def readSLHAFile(spcfilename, **kwargs):
    """
    Read an SLHA file, returning dictionaries of blocks and decays.

    Other keyword parameters are passed to readSLHA.
    """
    f = open(spcfilename, "r")
    rtn = readSLHA(f.read(), kwargs)
    f.close()
    return rtn


#def writeSLHAFile(spcfilename, blocks, decays, listBlocks, listDecays, **kwargs):
def writeSLHAFile(spcfilename, function, **kwargs):
    """
    Write an SLHA file from the supplied blocks and decays dicts.

    Other keyword parameters are passed to writeSLHA.
    """
    f = open(spcfilename, "w")
    f.write(writeSLHA(function, kwargs))
#    f.write(writeSLHA(blocks, decays, listBlocks, listDecays, kwargs))
    f.close()


def readISAWIGFile(isafilename, **kwargs):
    """
    Read a spectrum definition from a file in the ISAWIG format, returning
    dictionaries of blocks and decays. While this is not an SLHA format, it is
    informally supported as a useful mechanism for converting ISAWIG spectra to
    SLHA.

    Other keyword parameters are passed to readSLHA.
    """
    f = open(isafilename, "r")
    rtn = readISAWIG(f.read(), kwargs)
    f.close()
    return rtn


def writeISAWIGFile(isafilename, blocks, decays, **kwargs):
    """
    Write an ISAWIG file from the supplied blocks and decays dicts.

    Other keyword parameters are passed to writeISAWIG.
    """
    f = open(isafilename, "w")
    f.write(writeISAWIG(blocks, decays, kwargs))
    f.close()


def writeISAJETFile(isafilename, blocks, decays, **kwargs):
    """
    Write an ISAJET file from the supplied blocks and decays dicts (see writeISAJET).

    Other keyword parameters are passed to writeISAJET.
    """
    f = open(isafilename, "w")
    f.write(writeISAWIG(blocks, decays, kwargs))
    f.close()



###############################################################################
## Main function for module testing


if __name__ == "__main__":
    import sys
    for a in sys.argv[1:]:
        if a.endswith(".isa"):
            blocks, decays = readISAWIGFile(a)
        else:
            blocks, decays, listBlocks, listDecays = readSLHAFile(a)

        for bname, b in sorted(blocks.iteritems()):
            print b
            print

        print blocks.keys()

        print blocks["MASS"].entries[25]
        print

        for p in sorted(decays.values()):
            print p
            print

        print writeSLHA(blocks, decays, listBlocks, listDecays, ignorenobr=True)
