#! /usr/bin/env python

# derive exclusions for the tested scenario from the signal regions of an
# analysis implemented using MadAnalysis 5 (MA5), using a simplified procedure
# and the CLs prescription

# version 1.0 (July 9, 2014)
# made by Beranger Dumont
# based on toy MC code by Benjamin Fuks, Chris Wymant and Sam Bein

# takes as input:
# -- an XML analysis_name.info file that should be present in the same directory
# as the analysis code, i.e. Build/SampleAnalyzer/User/Analyzer/
# -- the SAF files for the cutflows of the tested scenario of interest in all
# signal regions, where the information on acceptance*efficiency can be found
# -- the cross section for the tested scenario (if not given as argument in
# command-line, it is looked for in the SAF file analysis_name.saf

# returns the output on the screen and print basic results
# in the file analysis_name_[run number].out in
# Output/benchmark_point/


# path to the directory where the analysis code and output is present
# (this is the directory created when running MA5 in expert mode,
# containing the "Build", "Input" and "Output" directories)
analysis_path = "./"

# The number of Poisson distributions we consider (each one effectively being
# a toy experiment with its own certain prediction for the background):
numtoyexperiments = 100000

#
# the user is not supposed to modify the code below this line
#

import os, sys
try:
    from lxml import ET
except:
    import xml.etree.ElementTree as ET
import scipy.stats

def usage():
    print 'Usage: ./exclusion_CLs.py analysis_name benchmark_point ' + \
      '[run_number] [cross section in pb]'
    print 'Example: ./exclusion_CLs.py cms_sus_13_011 T2tt_650_50.txt ' + \
      '0 0.014'
    print 'Default value of run number is 0'
    print 'If the cross section is not given, it is taken from the MA5 output'

def listdirectory(path):
    lfile=[]
    for root, dirs, files in os.walk(path): 
        for i in files: 
            lfile.append(os.path.join(root, i))
    return lfile

bestSR = ""

def clean_name(str):
    # based on CleanName() from
    # ./tools/SampleAnalyzer/Process/Core/SampleAnalyzer.cpp
    # in MA5 v1.1.11beta4
    str = str.replace("/",  "_slash_")
    str = str.replace("->", "_to_")
    str = str.replace(">=", "_greater_than_or_equal_to_")
    str = str.replace(">",  "_greater_than_")
    str = str.replace("<=", "_smaller_than_or_equal_to_")
    str = str.replace("<",  "_smaller_than_")
    str = str.replace(" ",  "_")
    str = str.replace(",",  "_")
    str = str.replace("+",  "_")
    str = str.replace("-",  "_")
    str = str.replace("(",  "_lp_")
    str = str.replace(")",  "_rp_")
    return str

def CLs(NumObserved, ExpectedBG, BGError, SigHypothesis, NumToyExperiments):
    # generate a set of expected-number-of-background-events, one for each toy
    # experiment, distributed according to a Gaussian with the specified mean
    # and uncertainty
    ExpectedBGs = scipy.stats.norm.rvs(loc=ExpectedBG, \
    scale=BGError, size=NumToyExperiments)

    # Ignore values in the tail of the Gaussian extending to negative numbers
    ExpectedBGs = [value for value in ExpectedBGs if value > 0]

    # For each toy experiment, get the actual number of background events by
    # taking one value from a Poisson distribution created using the expected
    # number of events.
    ToyBGs = scipy.stats.poisson.rvs(ExpectedBGs)
    ToyBGs = map(float, ToyBGs)

    # The probability for the background alone to fluctutate as LOW as
    # observed = the fraction of the toy experiments with backgrounds as low as
    # observed = p_b.
    # NB (1 - this p_b) corresponds to what is usually called p_b for CLs.
    p_b = scipy.stats.percentileofscore(ToyBGs, NumObserved, kind='weak')*.01

    # Toy MC for background+signal
    ExpectedBGandS = [expectedbg + SigHypothesis for expectedbg in ExpectedBGs]
    ToyBplusS = scipy.stats.poisson.rvs(ExpectedBGandS)
    ToyBplusS = map(float, ToyBplusS)

    # Calculate the fraction of these that are >= the number observed,
    # giving p_(S+B). Divide by (1 - p_b) a la the CLs prescription.
    p_SplusB = scipy.stats.percentileofscore(ToyBplusS, NumObserved, kind='weak')*.01
    
    return 1.-(p_SplusB / p_b) # 1 - CLs

def exclusion_check(crosssection):
    global bestSR

    if len(signalregions) > 1:
        # if more than one signal region, decide which signal regions
        # yields the best expected limit
        limit = -1.
        for SR in signalregions:
            nsignal = crosssection * lumi * 1000. * signalregions[SR]["acceff"]
            nb = signalregions[SR]["nb"]
            deltanb = signalregions[SR]["deltanb"]

            limitSR = CLs(nb, nb, deltanb, nsignal, numtoyexperiments)
            if limitSR > limit:
                bestSR = SR
                limit = limitSR
    else:
        bestSR = signalregions.keys()[0]

    nsignal = crosssection * lumi * 1000. * signalregions[bestSR]["acceff"]
    nobs = int(signalregions[bestSR]["nobs"])
    nb = signalregions[bestSR]["nb"]
    deltanb = signalregions[bestSR]["deltanb"]

    print 'The best expected signal region is "' + bestSR + '".'
    print 'It has: nobs = ' + str(nobs) + ', nb = ' + str(nb) + ' \pm ' + \
      str(deltanb) + ', nsignal = ' + str(round(nsignal,2)) + '.'

    return CLs(nobs, nb, deltanb, nsignal, numtoyexperiments)

def exclusion_check95(crosssection):
    return exclusion_check(crosssection)-0.95

if analysis_path[-1] == "/":
    analysis_path = analysis_path[:-1]

# at least one argument, check if asking for help
if len(sys.argv) < 3 or (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
    usage()
    sys.exit()

analysis_name = sys.argv[1]
bench_name = sys.argv[2]
if len(sys.argv) > 3:
    run_number = sys.argv[3]
else:
    run_number = "0"
if len(sys.argv) > 4:
    try:
        xsection = float(sys.argv[4])
    except ValueError:
        print 'Invalid cross section given as command-line argument: "' + \
          sys.argv[4]+ '".'
        sys.exit()
else:
    xsection = 0



######################################
# first, read the XML .info file
# and fill the variable lumi
# and the dictionary signalregions
######################################

lumi = 0 # integrated luminosity, in fb^-1
signalregions = {}

analysisinfo_path = analysis_path + "/Build/SampleAnalyzer/User/Analyzer/" + \
    analysis_name + ".info"

try:
    info_input = open(analysisinfo_path)
except IOError as e:
    print 'I/O error({0}): {1}'.format(e.errno, e.strerror)
    print 'Cannot open the XML info file "' + analysisinfo_path + '".'
    sys.exit()

info_tree = ET.parse(info_input)
info_input.close()

info_root = info_tree.getroot()
if info_root.tag != "analysis":
    print 'Invalid XML info file "' + analysisinfo_path + '".'
    print 'Root tag should be <analysis>, not <' + info_root.tag + '>.'
    sys.exit()
if info_root.attrib["id"] != analysis_name:
    print 'Invalid XML info file "' + analysisinfo_path + '".'
    print 'Analysis id in root tag <analysis> should be "' + analysis_name + \
          '", not "' + info_root.attrib["id"] + '".'
    sys.exit()


for child in info_root:
    # for <lumi> tag
    if child.tag == "lumi":
        if lumi != 0:
            print 'Warning: redefinition of the luminosity in the ' + \
              'XML info file "' + analysisinfo_path + '".'

        try:
            lumi = float(child.text)
        except TypeError: # empty tag is of type NULL
            lumi = 0
        except ValueError:
            print 'Invalid XML info file "' + analysisinfo_path + '".'
            print 'The value of the <lumi> tag is not a number.'
            sys.exit()

    # for the <region> tags
    # if no type is specified, assumed to be a signal region
    if child.tag == "region" and \
      ("type" not in child.attrib or child.attrib["type"] == "signal"):
        if "id" not in child.attrib:
            print 'Invalid XML info file "' + analysisinfo_path + '".'
            print 'Presence of <region> tags with no id attribute.'
            sys.exit()

        regionid = child.attrib["id"]

        if regionid in signalregions:
            # a <region> tag with the same id has already been defined
            print 'Invalid XML info file "' + analysisinfo_path + '".'
            print 'A region with id="' + regionid + ' is defined ' + \
              'multiple times.'
            sys.exit()

        signalregions[regionid] = {"acceff": 0} # initialize efficiency to 0
        for rchild in child:
            if rchild.tag in ["nobs", "nb", "deltanb"]:
                ntag = rchild.tag
                if ntag in signalregions[regionid]:
                    print 'Warning: redefinition of <' + ntag + '> in the ' + \
                      'region "' + \
                      regionid + '" of the XML info file "' + \
                      analysisinfo_path + '".'

                try:
                    signalregions[regionid][ntag] = float(rchild.text)
                except TypeError: # empty tag is of type NULL
                    signalregions[regionid][ntag] = 0
                except ValueError:
                    print 'Invalid XML info file "' + analysisinfo_path + '".'
                    print 'The value of the <' + ntag + '> tag in region "' + \
                     regionid + '" is not a number.'
                    sys.exit()

######################################
# then, read the SAF files
# for the cutflows
# generated by MA5
######################################

analysisinfo_path = analysis_path + "/Output/" + \
    bench_name + "/" + analysis_name + "_" + run_number + "/Cutflows"

listdir = listdirectory(analysisinfo_path)

if not listdir:
    print 'The directory "' + analysisinfo_path + '" containing the SAF' + \
      ' files for the cutflows cannot be listed.'
    sys.exit()

for file in listdir:
    # signal region (as defined in the XML info file)
    # associated with the cutflow file
    assoc_SR = ""

    if file[-4:] != ".saf":
        continue
    SRname = file.split("/")[-1][:-4]

    for regionid in signalregions:
        for sr in regionid.split(";"):
            if clean_name(sr) == SRname:
                assoc_SR = sr

    # if there is no signal region found associated with the SAF file
    if assoc_SR == "":
        print 'Warning: no region found associated with the SAF file "' + \
          file + '"; will be skipped.'
        continue

    # otherwise, read the acceptance times efficiency from the SAF file
    try:
        SAF_cutflow = open(file)
    except IOError as e:
        print 'I/O error({0}): {1}'.format(e.errno, e.strerror)
        print 'Cannot open the XML info file "' + file + '".'
        sys.exit()

    in_initialcounter = False
    in_counter = False
    initialnum = 0
    finalnum = 0
    for line in SAF_cutflow:
        line = line.rstrip("\n").strip()
        if line[:16] == "<InitialCounter>":
            in_initialcounter = True
            continue
        if line[:17] == "</InitialCounter>":
            in_initialcounter = False
            continue
        if line[:9] == "<Counter>":
            in_counter = True
            continue
        if line[:10] == "</Counter>":
            in_counter = False
            continue

        if in_initialcounter and line[-14:] == "sum of weights":
            try:
                initialnum = float(line.split()[0])
            except:
                print 'Invalid SAF file "' + file + '".'
                print 'The initial number of events cannot be read.'
                sys.exit()

        if in_counter and line[-14:] == "sum of weights":
            try:
                finalnum = float(line.split()[0])
            except:
                print 'Invalid SAF file "' + file + '".'
                print 'The number of events cannot be read.'
                sys.exit()

    if initialnum == 0:
        print 'Invalid SAF file "' + file + '".'
        print 'The number of events cannot be read.'
        sys.exit()

    for regionid in signalregions:
        # for each SR id as defined in the XML info file...
        for individualSR in regionid.split(";"):
            # ...look for each individual SR...
            if individualSR == assoc_SR:
                # ... if it matches with the SAF file read, add the
                # acceptance*efficiency to existing value
                signalregions[regionid]["acceff"] += finalnum / initialnum

######################################
# then, get the cross section info
# if not given as command-line
# argument, look into SAF file
######################################

if xsection == 0: # not given as command-line argument
    mainSAF_path = analysis_path + "/Output/" + \
    bench_name + "/" + bench_name + ".saf"

    try:
        mainSAF = open(mainSAF_path)
    except IOError as e:
        print 'I/O error({0}): {1}'.format(e.errno, e.strerror)
        print 'Cannot open the XML info file "' + mainSAF_path + '".'
        sys.exit()

    nline = 0
    for line in mainSAF:
        line = line.rstrip("\n").strip()
        if line[:18] == "<SampleGlobalInfo>":
            nline = 1
            continue
        if nline > 0:
            nline += 1
        if nline == 3:
            try:
                xsection = float(line.split()[0])
            except:
                print 'Invalid SAF file "' + mainSAF_path + '".'
                print 'The cross section cannot be read.'
                sys.exit()

            if xsection <= 0.:
                print 'Invalid cross section of ' + str(xsection) + ' pb.'
                print 'The cross section cannot be read from the SAF file "' + \
                  mainSAF_path + '".'
                sys.exit()
            break

######################################
# now, we have the information on the
# cross section (variable xsection),
# the luminosity (variable lumi),
# the number of background events,
# observed events, and the acceptance
# times efficiency (in the dictionary
# signalregions)
#
# we can proceed with the exclusion
######################################

if xsection > 0:
    # first, decide which signal regions yields the best expected limit
    final_limit = exclusion_check(xsection)

    print '\nThis signal is excluded at the ' + str(round(final_limit*100.,1)) + \
      '% CL (CLs=' + str(round(1-final_limit,3)) + ').'
else:
    # if a  negative cross section is given as input,
    # the code is looking for the cross section that is excluded at 95% CL
    # using a root-finding algorithm

    print 'Negative cross section is given.'
    print 'Will look for the cross section that is excluded at 95% CL'
    print 'using a root-finding algorithm.\n'

    # need to tune the lower and upper bound (corresponding to a cross section
    # that we know is not excluded or excluded, respectively)
    lowerb = 1. # 1 pb
    upperb = 1. # 1 pb
    while exclusion_check95(lowerb) > 0.:
        lowerb *= 0.1
    while exclusion_check95(upperb) < 0.:
        upperb *= 10.

    print '\nlower and upper bounds for the root-finding algorithm' + \
      ' have been found: [%.2e %.2e] pb\n' % (lowerb, upperb)

    final_limit = scipy.optimize.brentq(exclusion_check95, lowerb, upperb)

    print '\nThe excluded cross section at 2 sigma is %.5E pb.' % final_limit

# finally, write the results in an output file
output_path = analysis_path + "/Output/" + \
    bench_name + "/" + analysis_name + "_" + run_number + ".out"

try:
    output = open(output_path, "w")
except IOError as e:
    print 'I/O error({0}): {1}'.format(e.errno, e.strerror)
    print 'Cannot create the output file "' + output_path + '".'
    sys.exit()

output.write(bestSR+'\n'+str(final_limit))
output.close()
