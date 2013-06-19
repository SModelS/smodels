#!/usr/bin/python
import set_path, argparse, types, os

argparser = argparse.ArgumentParser(description='simple tool to generate a latex table with all analysis used') 
argparser.add_argument ( '-o', '--output', nargs='?', help='output file', type=types.StringType, default='tab.tex')
argparser.add_argument ( '-D', '--Dir', nargs='?', help='figure folder', type=types.StringType, default='figures/')
args=argparser.parse_args()

from Experiment import SMSAnalysisList, SMSAnalysisFactory
from Tools import TexTable

DoFactory = True
#Creat analyses list:
if DoFactory:
  ListOfAnalyses = SMSAnalysisFactory.load()
else:  
  ListOfAnalyses = SMSAnalysisList.load()

if not os.path.isdir(args.Dir):
  os.system("mkdir "+args.Dir)

#Generate table:
TexTable.GenAnalysisTable(ListOfAnalyses, texfile=args.output, fig_dir=args.Dir)

