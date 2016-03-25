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
    parser.add_argument(
        "--dryrun",   
        help='print out a list of changes which would be made but don\'t change any files ',action="store_true")

    args = parser.parse_args()

    print "headerscan: Making a list of local header files."
    print "===============================================\n"

    localheaders = glob.glob('*.h')
    
    pprint(localheaders)

    print "===============================================\n"

    print "headerscan: Scanning over files in ./src.\n"

    p = re.compile('#include<(?P<headerfile>\S+)>(?P<endmatter>.*)')
    
    for line in fileinput.input(glob.glob("pdcode.c"),inPlace=True,backup='.backup'):

        m = p.match(line)
        if m:
            
            if (m.group('headerfile') in localheaders):

                lst = ['#include"',m.group('headerfile'),'"',m.group('endmatter')]                  if args.dryrun:
                    print "would replace ", line.rstrip(), "with ","".join(lst)
                    sys.stdout.write(line)
                else
                    sys.stdout.write("".join(lst))

        else:

            sys.stdout.write(line)

    
