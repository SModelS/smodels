# the pmssm walk ding
# (c) Wolfgang Waltenberger, 2012

import sys, os, subprocess

progname=sys.argv[0]
test=""
if progname.find("pmssmtest")>0:
  test="test"
basedir="/afs/hephy.at/user/w/walten/pmssm"+test
wwwdir="/afs/hephy.at/user/w/walten/www/pmssm"+test
url="http://www.hephy.at/user/walten/pmssm"+test
cgi="http://www.hephy.at/user/walten/cgi-bin"
uploaddir="%s/upload" % basedir
datadir="%s/data" % basedir
logdir="%s/log" % basedir
slhadir="%s/slha" % wwwdir
plotsdir="%s/plots" % wwwdir
feyndir="%s/feyn" % wwwdir
templatedir="%s/templates" % basedir
installdir="%s/install" % basedir
etcdir="%s/etc" % basedir
rootdir="%s/install/lib/root" % basedir
forceData=False ## force the production of data files
forceHtml=False ## force the creation of html files
pNMSSM=False ## is this the pNSSM case?
printmass=False ## show the masses in the rulerplot
darkmatter=False ## compute dark matter quantities

verbose=False ## show the logs
noCache=False ## use the cached plots and everything
description=False ## add description to pmssm parameter form
lhe1000=False ## generate 1000 events
considerSezen=False ## consider even looking for sezen info?
sezen=False ## look for sezen info?
block=False ## block CMS-specific info?
lhe10000=False ## generate 10000 events

ModelTag="pMSSM"
if pNMSSM: ModelTag="pNMSSM"

similarSMSes= { "TGQ3": "TGQ", "T5zzoff": "T5zz",
  "TGQ3inter": "TGQ", "TChiwzincl": "TChiwz",
  "T6ttwzoff": "T6ttwz", "T6ttzzoff": "T6ttzz",
  "TChiwwincl": "TChiww", "T6bbzzoff": "T6bbzz",
  "T6ttwwoff": "T6ttww", "T6ttzzoff": "T6ttzz",
  "T2ttincl": "T2tt", "T2bbincl": "T2bb",
  "TChizzincl": "TChizz", "TGQinter":"TGQ", "T5zw":"T5wz",
  "TChiC1C1": "TChiww", "TChiC1C1+": "TChiww",
  "TChiN2C1": "TChiwz", "TChiN2N2": "Tchizz",
  "T5CC": "T5ww", "T5NN": "T5zz", "T5CN": "T5wz",
  "T2general": "T2", "TChiwh+": "TChiwh",
  "T6bbWW5GeV": "T6bbWW", "T6bbWW20GeV": "T6bbWW",
  "TChiN2C1?": "TChiwz", 
  "TChiN2C1+": "TChiwz", 
  "T6bb*": "T2bb", 
  "TChiN2C1": "TChiwz", 
  "TChiN2N2": "TChizz", 
  "TChiC1C1": "TChiww", 
  "TChiC1C1?": "TChiww", 
  "TChiC1N1+": "TChiC1N1", 
  "T4bbw": "T2bb", 
  "T4bbz": "T2bb", 
  "T4ttN": "T2tt", 
  "T6ttNN": "T2tt", 
  "T6ttWbWb": "T2tt", 
  "T2ttWbWb": "T2tt", 
  "TChiC1N1": "TChiwn",
  "TChiC1C1": "TChiww",
  "TChiC2C2": "TChiww",
  "TChiN2N1*": "TChizn",
  "TSlepSnu+": "TSlepSnu",
  "TSlepSlep+": "TSlepSlep",
  "TSnuSnu+": "TSnuSnu",
  "TChiSlepSlep+": "TChiSlepSlep",
}

smsDescriptions= {
  "T6bbWW5GeV": "Same as T6bbWW but with chargino-LSP mass splitting of 5 GeV (see ATLAS-CONF-2013-001)",
  "T6bbWW20GeV": "Same as T6bbWW but with chargino-LSP mass splitting of 20 GeV (see ATLAS-CONF-2013-001)",
  "TChiwh": "C1 -> W LSP, N2 -> h LSP",
  "TChiwh+": "production of C1 and N2, multiple higgses and W bosons",
}

def log ( text, level="log" ):
  printit=True
  if not verbose and level=="log":
    printit=False
  if level!="log" and level!="warn" and level!="info":
    print "Unknown log level",level,"<br>"
  color,endcolor="",""
  if level=="info":
    color="<font color=darkblue>"
    endcolor="</font>"
  if printit: print "%s<%s>[%s] %s</%s>%s<br>" % ( color, level, level, text, level, endcolor )

def fixme ( text ):
  print "<fixme>FIXME: %s</fixme><br>" % text

parameterdescription= {
  ### these are the names and the explanations for the pmssm parameters
  "sqrts": "center-of-mass energy",
  "tanbeta": "tan_beta(m_Z)",
  "mtop": "m_top(pole)",
  "alphas": "alpha_s(MZ)",
  "M1": "Gaugino mass M_1",
  "M2": "Gaugino mass M_2",
  "M3": "Gaugino mass M_3",
  "Atop": "Trilinear coupling A_top",
  "Abottom": "Trilinear coupling A_bottom",
  "Atau": "Trilinear coupling A_tau",
  "Higgsino": "Higgsino mass mu",
  "mA0": "Pseudoscalar Higgs pole mass m_A0",
  "mslepl": 'Left 1st gen. scalar lepton mass',
  "mstaul": 'Left 3rd gen. scalar lepton mass',
  "mslepr": 'Right scalar electron mass',
  "mstaur": "Right scalar tau mass",
  "msquarkl": "Left 1st gen. scalar quark mass",
  "mstopl": "Left 3rd gen. scalar quark mass",
  "msupr": "Right scalar up mass",
  "mstopr": "Right scalar top mass",
  "msdownr": "Right scalar down mass",
  "msbottomr": "Right scalar bottom mass",
}

if pNMSSM:
  parameterdescription["lambda"]="?"
  parameterdescription["kappa"]="?"
  parameterdescription["alambda"]="?"
  parameterdescription["akappa"]="?"


parameterhtml= {
  ### these are the names and the explanations for the pmssm parameters
  #"sqrts": "&radic;s",
  "sqrts": '&radic;<span style="text-decoration:overline;">s</span>',
  "tanbeta": "tan &beta;",
  "mtop": "m(t)",
  "alphas": "&alpha;<sub>s</sub>",
  "M1": "M<sub>1</sub>",
  "M2": "M<sub>2</sub>",
  "M3": "M<sub>3</sub>",
  "Atop": "A(t)",
  "Abottom": "A(b)",
  "Atau": "A(&tau;)",
  "Higgsino": "&mu;",
  "mA0": "m(A<sub>0</sub>)",
  "mslepl": 'm(&tilde;l<sub>L</sub>)',
  "mstaul": 'm(&tilde;&tau;<sub>L</sub>)',
  "mslepr": 'm(&tilde;l<sub>R</sub>)',
  "mstaur": "m(&tilde;&tau;<sub>R</sub>)",
  "msquarkl": "m(&tilde;q<sub>L</sub>)",
  "mstopl": "m(&tilde;t<sub>L</sub>)",
  "msupr": "m(&tilde;u<sub>R</sub>)",
  "mstopr": "m(&tilde;t<sub>R</sub>)",
  "msdownr": "m(&tilde;d<sub>R</sub>)",
  "msbottomr": "m(&tilde;b<sub>R</sub>)",
  "lambda": "&lambda;",
  "kappa": "&kappa;",
  "alambda": "a&lambda;",
  "akappa": "a&kappa;",
}

parameterlatex= {
  ### these are the names and the explanations for the pmssm parameters
  "tanbeta": "\\tan \\beta",
  "mtop": "m(t)",
  "alphas": "\\alpha_{s}",
  "sqrts": "\\sqrt{s}",
  "M1": "M_{1}",
  "M2": "M_{2}",
  "M3": "M_{3}",
  "Atop": "A(t)",
  "Abottom": "A(b)",
  "Atau": "A(\\tau)",
  "Higgsino": "\\mu",
  "mA0": "m(A_0)",
  "mslepl": 'm(\\tilde{l}_{L})',
  "mstaul": 'm(\\tilde{\\tau}_{L})',
  "mslepr": 'm(\\tilde{l}_{R})',
  "mstaur": "m(\\tilde{\\tau}_{R})",
  "msquarkl": "m(\\tilde{q}_{L})",
  "mstopl": "m(\\tilde{t}_{L})",
  "msupr": "m(\\tilde{u}_{R})",
  "mstopr": "m(\\tilde{t}_{R})",
  "msdownr": "m(\\tilde{d}_{R})",
  "msbottomr": "m(\\tilde{b}_{R})",
  "lambda": "\\lambda",
  "kappa": "\\kappa",
  "alambda": "a\\lambda",
  "akappa": "a\\kappa",
}


parameternumbers= {
  ### these are the numbers of the parameters in the MINPAR block,
  ### plus 1000 to distinguish from EXTPAR
  1003: "tanbeta",
  ### these are the numbers of the parameters in the SMINPUTS block,
  ### plus 2000 to distinguish from EXTPAR
  2006: "mtop",
  2003: "alphas",
  9999: "sqrts",
  ### these are the numbers of the parameters in the EXTPAR block
  1:"M1",
  2:"M2",
  3:"M3",
  11:"Atop",
  12:"Abottom",
  13:"Atau",
  23:"Higgsino",
  26:"mA0",
  31:"mslepl",
  32:"mslepl",
  33:"mstaul",
  34:"mslepr",
  35:"mslepr",
  36:"mstaur",
  41:"msquarkl",
  42:"msquarkl",
  43:"mstopl",
  44:"msupr",
  45:"msupr",
  46:"mstopr",
  47:"msdownr",
  48:"msdownr",
  49:"msbottomr",
  61:"lambda",
  62:"kappa",
  63:"alambda",
  64:"akappa",
  65:"mueff"
}

floatparams=[ "alphas", "mtop", "lambda", "kappa" ] ## , "alambda", "akappa" ]

parameterdefaults= {
  ### some sensible defaults
  "sqrts": 7,
  "tanbeta": 30,
  "alphas": .1184,
  "mtop": 172.9,
  "M1": 300,
  "M2": 400,
  "M3": 1400,
  "Atop": 800,
  "Abottom": 0,
  "Atau": 0,
  "Higgsino": 400,
  "mA0": 1000,
  "mslepl": 600,
  "mstaul": 1700,
  "mslepr": 1500,
  "mstaur": 1000,
  "msquarkl": 1250,
  "mstopl": 1750,
  "msupr": 1500,
  "mstopr": 1750,
  "msdownr": 1250,
  "msbottomr": 1000,
  "lambda": .7,
  "kappa": .05,
  "alambda": 1280,
  "akappa": 0
}

def num_in_base(val, base=62, min_digits=1, complement=False,
       digits="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"):
    """Convert number to string in specified base
       If minimum number of digits is specified, pads result to at least
       that length.
       If complement is True, prints negative numbers in complement
       format based on the specified number of digits.
       Non-standard digits can be used. This can also allow bases greater
       than 36.
    """
    if base < 2: raise ValueError("Minimum base is 2")
    if base > len(digits): raise ValueError("Not enough digits for base")
    # Deal with negative numbers
    negative = val < 0
    val = abs(val)
    if complement:
        sign = ""
        max = base**min_digits
        if (val > max) or (not negative and val == max):
            raise ValueError("Value out of range for complemented format")
        if negative:
            val = (max - val)
    else:
        sign = "-" * negative
    # Calculate digits
    val_digits = []
    while val:
        val, digit = divmod(val, base)
        val_digits.append(digits[digit])
    result = "".join(reversed(val_digits))
    leading_digits = (digits[0] * (min_digits - len(result)))
    return sign + leading_digits + result

def createKey ( dic ):
  """ create a unique key for a specific choice of pmssm parameters. """
  haskeys={}
  keys=parameterdescription.keys()
  keys.sort()
  hasnokeysatall=True
  for i in keys:
    haskeys[i]=False
  masterkey=""
  for (param,value) in dic.items():
    if param=="comment": continue
    if param=="darkmatter": continue
    if not haskeys.has_key ( param ):
      print "[%s] unknown key %s" % ( ModelTag, param )
      sys.exit(0)
    haskeys[param]=True
    hasnokeysatall=False

  for (key,value) in haskeys.items():
    if not value:
      Def=parameterdefaults[key]
      if not hasnokeysatall:
        ## issue the warning only if *some* values are missing!
        log("did not supply value for: %s. Will use default of %s" \
              % ( key, Def ), "warn" )
      dic[key]=Def

  for param in keys:
    value=dic[param]
    ##indx=keys.index ( param )
    if not param in floatparams:
    #if param!="mtop" and param!="alphas":
      value=int(round(value))
    if value<0: value=100000+abs(value)
    if ( int(float(value))== value ):
      masterkey += str (int(float(value)))
    else:
      masterkey += str (int(10000*float(value)))
    if masterkey.find(".")>-1:
      print "[PMSSM.py] Error: found ``.'' in masterkey.<br>"
  ret=num_in_base ( long(masterkey) )
  return ret

def produceLesHouchesInFile ( key, params, outfile ):
  """ this file takes a masterkey, a params dictionary, and the
    name of the output file and produces a bare slha file,
    containing the input parameters and the SM parameters,
    but not the particles masses and decays. """
  try:
    f=open("%s/template.%s.in" % ( templatedir, ModelTag ) )
    lines=f.readlines()
    #filename="LH_%s.in" % (key)
    #outfile="%s/%s" % (datadir,filename)
    w=open( outfile,"write")
    for line in lines:
      for (key,value) in params.items():
        line=line.replace ( "PMSSM_%s" % key, str(value) )
      w.write ( line )
    w.close()
    log ( "produced %s" % outfile )
  except Exception,e:
    ##e="unknown"
    log ( "failed: %s" % e )

def runSPheno ( key ):
  """ run spheno """
  import popen2
  # lhfile="%s/LH_%s.in" % ( datadir, key )
  ## slhafile="%s/%s.slha" % (datadir, key )
  lhfile="../data/LH_%s.in" % key
  slhafile="../data/%s.slha" % key
  slhafileabs="%s/%s.slha" % (datadir, key )
  cmd="cd %s; %s/SPheno %s %s; cd %s" % \
    (datadir, installdir,lhfile,slhafile, installdir)
  ##cmd="./SPheno %s" % (lhfile)
  ## cmd="./SPheno"
  log ( "now run %s" % cmd )
  pp=popen2.popen3( cmd )
  lines=pp[2].readlines()
  for line in lines: log ( "[SPheno]"+line  )
  # cmd="cp %s %s/%s.slha" % ( lhfile, slhadir, key )
  cmd="cp %s %s" % ( slhafileabs, slhadir )
  pp1=popen2.popen3 ( cmd )
  lines=pp1[2].readlines()
  for l in lines: log ( "[cp to data] %s<br>" % l )
  #cmd="cp %s %s" % ( lhfile, slhafile )
  #log ( cmd )
  #popen2.popen3 ( cmd )
  cmd="cp %s %s/%s.slha" % ( slhafileabs, slhadir, key )
  log ( cmd )
  pp2=popen2.popen3 ( cmd )
  log ( "done, new slha file: %s" %  slhafile )
  #log ( "pp2="+str(pp2) )
  lines=pp2[2].readlines()
  for l in lines: log ( "[cp to www] %s<br>" % l )
  f=open("%s/Messages.out" % datadir )
  lines=f.readlines()
  f.close()
  for line in lines:
    if line.find("CP violation"): continue
    log ( line, "warn" )

def produceParameterFile ( key, picklefile, params ):
  """ dump the parameter dictionary in a text file.
      return value: has the pickle file already existed?
  """
  existed=os.path.exists ( picklefile )
  f=open( picklefile,"w")
  f.write ( str(params)+"\n" )
  # pickle.dump ( params, f )
  f.close()

def runNMSSM ( key ):
  """ run the NMSSM tools """
  import runNMSSMTools
  runNMSSMTools.produce ( key, verbose  )

def runPythiaLHE ( n, slhafile, masterkey, logMe=False, sqrts=7 ):
  """ run pythia_lhe with n events, at sqrt(s)=sqrts.
      slhafile is inputfile
      datadir is where this all should run,
      installdir is where pythia_lhe is to be found.  """
  import commands
  log ( "try to run pythia_lhe at sqrt(s)=%d with %d events" % (sqrts,n) )
  o=commands.getoutput ( "cp %s %s/fort.61" % ( slhafile, datadir ) )
  if len(o)>0:
    print "[PMSSM.py] runPythiaLHE error",o
  f=open(etcdir+"/external_lhe.template") 
  lines=f.readlines()
  f.close()
  g=open(datadir+"/external_lhe.dat","write")
  for line in lines:
    out=line.replace("NEVENTS",str(n)).replace("SQRTS",str(1000*sqrts))
    g.write ( out )
  g.close()
  o=commands.getoutput ( "cd %s; %s/pythia_lhe < %s/external_lhe.dat" % \
     ( datadir, installdir, datadir ) )
  lines=o.split( "\n" )
  xsecfb=None
  flog=None
  if logMe:
    flog=open("%s/pythia_%s.log" % ( logdir, masterkey ), "write" )
  for line in lines:
    if logMe: flog.write ( line+ "\n" )
    if line.find("All included subprocesses")>0:
      try:
        xsecfb=float(line[67:78].replace("D","E"))*10**12
      except Exception,e:
        print "[ResultsTables.py] Exception",e,"xsecfb=",line[67:78]
        print "  `-- line=",line
        print "  `-- masterkey=",masterkey
  if logMe: flog.close()
  return { "xsecfb": xsecfb }
