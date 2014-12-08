#!/usr/local/bin/ipython

import libpl.pdstor as pdstor
import libpl.pdcode as pdcode
import spherogram as sgm
import sys, getopt 
from libpl.util.edit import pd_edit

def main(argv):

    outputfile = ''
    try:
      opts, args = getopt.getopt(argv,"h:o:",["ofile="])
    except getopt.GetoptError:
      print 'sketchpd.py -o <outputfile>'
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
        print 'sketchpd.py -o <outputfile>'
        sys.exit()
      elif opt in ("-o", "--ofile"):
        outputfile = arg
    
    print 'sketchpd started. Output file is "',outputfile
    new_code = pd_edit()
    with open(outputfile,'w') as f:
        new_code.write(f)
        f.close()

if __name__ == "__main__":
   main(sys.argv[1:])
   
    
    
