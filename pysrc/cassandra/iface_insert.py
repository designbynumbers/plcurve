#!/usr/bin/env python

""" An example python script which uses the ColdInterface class to load in some pdcodes
"""
from iface import ColdInterface

HOSTNAME = "localhost"
PDSTOR_DIR = "../../data/pdstors"
CROSSINGS = [3,4,5]
if __name__ == "__main__":
    cold = ColdInterface(HOSTNAME, dirloc=PDSTOR_DIR)
    print "Interface set up, inserting diagrams..."
    cold.insert_all_and_wait(CROSSINGS)
    print "Finished!"
