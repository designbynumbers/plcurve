#!/usr/bin/python
"""monogon_counter.py

author: Harrison Chapman
Counts monogon shadows in a pdstor file.
"""

from libpl.pdcode import PlanarDiagram

def monogon_in_diagram(pdc):
    for face in pdc.faces:
        if len(face) == 1:
            return True
    return False

def count_monogons(pdc):
    return sum(1 if len(face)==1 else 0 for face in pdc.faces)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Counts monogon shadows in a pdstor file."
    )
    parser.add_argument('--no-header', action="store_true",
                        help="do not try to read a header from the pdstor file")
    parser.add_argument('--no-links', action="store_true",
                        help="only count diagrams with a single component")
    parser.add_argument('files', metavar='file', nargs='+', type=argparse.FileType('r'),
                        help="files in which to count monogon shadows")

    args = parser.parse_args()

    for f in args.files:
        monogon_shadow_count = 0
        for pdc in PlanarDiagram.read_all(f, read_header=not args.no_header):
            if args.no_links and pdc.ncomps > 1:
                continue
            #monogon_shadow_count += count_monogons(pdc)
            if monogon_in_diagram(pdc):
                monogon_shadow_count += 1

        print("There are {} monogon shadows in the file {}".format(monogon_shadow_count, f.name))
