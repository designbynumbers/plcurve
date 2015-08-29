from libpl.pdcode import *
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file",help="input pdstor file")
parser.add_argument('-o','--outfile',help='output file',dest='ofile')
args = parser.parse_args()

if args.ofile != None:
    of = open(args.ofile,'w')
else:
    of = sys.stdout
      
with open(args.file) as infile:
    allpds = PlanarDiagram.read_all(infile,read_header=True)
    for pd in allpds:
        ed = pd.as_spherogram().view()
        of.write('{hash} {uid}\n'.format(hash=str(pd.hash),uid=str(pd.uid)))
        for arrow in ed.Arrows:
            of.write('\t {ar}\n'.format(ar=str(arrow)))
        of.write('\n')




