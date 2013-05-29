""" totally trivial stuff that is nice for testing """

green='\033[0;32m'
red='\033[0;31m'
reset='\033[;0m'

def ok ( A,B,verbose=True ):
  if A==B: return "%sok.%s" % ( green, reset )
  if verbose:
    return "%sfailed. [%s]%s" % ( red, B, reset )
  return "%sfailed. %s" % ( red, reset )

def convertROOTpdf(pdfin,pdfout,pngout):
  ''' take ROOT-pdf pdfin, create pdfout, pngout 
      NOT WORKING LIKE IT SHOULD'''
  from pdfrw import PdfReader, PdfWriter
  import os
  reader = PdfReader(pdfin)
  page, = reader.pages # read first page
  writer = PdfWriter()
  writer.addpage(page) ## adjust(page))
  writer.write(pdfout)
  os.system("convert %s %s" % (pdfout,pngout))
  return

