#!/usr/bin/python
import argparse
import tempfile
import os
import subprocess
import re
import sys
import atexit
import time
import glob
import shutil
import fileinput
from   pprint import pprint

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='scan for header files included with <header.h> instead of \"header.h\" ')
    parser.add_argument("-d","--dryrun",help="print out a list of changes which would be made but don\'t change any files",action="store_true")

    args = parser.parse_args()

    print "headerscan: Making a list of local header files."
    print "===============================================\n"

    localheaders = glob.glob('*.h')
    
    pprint(localheaders)

    print "===============================================\n"

    print "headerscan: Scanning over files in ./src.\n"

    p = re.compile('\s*#include\s*<(?P<headerfile>\S+)>(?P<endmatter>.*)')
    oldsrcfile = None
    
    for line in fileinput.input(glob.glob("*.c"),
                                inplace=(args.dryrun != None),backup='.backup'):

        if fileinput.filename() != oldsrcfile:
            oldsrcfile = fileinput.filename()
            sys.stderr.write("".join(["In file ",fileinput.filename()]))
                             
        m = p.match(line)
        if m:
            
            if (m.group('headerfile') in localheaders):

                lst = ['#include"',m.group('headerfile'),'"',m.group('endmatter'),"\n"]

                if args.dryrun:
                    sys.stderr.write("".join(["\twould replace ", line.rstrip(), " with ","".join(lst),"\n"]))
                    sys.stdout.write(line)
                else:
                    sys.stderr.write("".join(["\treplaced ",line.rstrip()," with ","".join(lst),"\n"]))
                    sys.stdout.write("".join(lst))

            else:
                sys.stdout.write(line)

        else:
            sys.stdout.write(line)

    
    print "\n\nheaderscan: diffing everything to show you the changes"
    print "------------------------------------------------------"

    for srcfile in glob.glob("*.c"):

        print("In file ",srcfile," ")
        subprocess.call(['diff',srcfile,"".join([srcfile,".backup"])])

     
