#!/usr/bin/env python
# check_pdstors.py
"""
Check that the PDStor files created by graph expansion are
equivalent to those which are created by CS of prime diagrams (duals
to simple quadrangulations)
"""

import os, os.path
RD_SRC_DIR = os.path.join(os.path.expanduser("~"),
                          "src", "randomdiagram", "src")
PD_DAT_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                          "..", "..", "data", "pdstors")

from libpl.pdstor import *

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        "Compare connected sum pdstors to graph expansion pdstors")
    parser.add_argument('ncross', type=int,
                        help="number of crossings")
    parser.add_argument('ncs', type=int, nargs="?", default=None,
                        help="number of connect sums (optional)")
    args = parser.parse_args()

    if args.ncs is not None:
        n, m = (args.ncross, args.ncs)
        cs_fname = "%scs%s.pdstor"%(n, m)
        ex_fname = "%s.pdstor"%(n + m,)
    else:
        n, m = (args.ncross, None)
        cs_fname = "%sb.pdstor"%(n,)
        ex_fname = "%s.pdstor"%(n,)
        
    with open(os.path.join(RD_SRC_DIR, cs_fname), "rb") as cs_f, \
         open(os.path.join(PD_DAT_DIR, ex_fname), "rb") as ex_f:
        print "Reading in connect sum file..."
        cs_pds = PDStorage.read(cs_f, read_header=True)
        print "done."
        print "Reading in graph expand file..."
        ex_pds = PDStorage.read(ex_f, read_header=True)
        print "done."

        if m is None:
            print "Checking that cs_pds == ex_pds..."
            equal = cs_pds == ex_pds
            print equal
            if not equal:
                "Error.. checking if cs_pds \subset ex_pds..."
                smaller = cs_pds < ex_pds
                print smaller
                "Error.. checking if cs_pds \supset ex_pds..."
                bigger = cs_pds > ex_pds
                print bigger
        else:
            print "Checking that cs_pds \subseteq ex_pds...",
            smallereq = cs_pds <= ex_pds
            print smallereq
            if not smallereq:
                "Error.. checking if cs_pds \supset ex_pds...",
                bigger = cs_pds > ex_pds
                print bigger

        

