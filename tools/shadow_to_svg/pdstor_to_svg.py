from shadow_to_svg import pd_to_svg
from libpl.pdcode import *
import os.path

def pdstor_to_svg(pdstor_f, read_header, subdirize=True, prefix="pd"):
    if subdirize:
        # Parse the header to get the count of shadows
        header = pdstor_f.readline()
        head_data  = pdstor_f.readline()
        _, nelts, _, _, nhashes = head_data.strip().split(" ")
        nclaimed, n_pd = (int(i) for i in nelts.split("/"))
        pdstor_f.seek(0)
        ndigits = len(str(n_pd))
        hunds = -1

    for i, pd in enumerate(PlanarDiagram.read_all(pdstor_f, read_header=read_header)):
        if subdirize:
            new_hunds = i/100
            si = str(i)
            si = "0"*(ndigits-len(si))+si
            subdir = os.path.join(prefix,*[si[:z+1] for z in xrange(len(si)-2)])
            if hunds != new_hunds:
                hunds = new_hunds
                os.makedirs(subdir)
            pd_to_svg(pd, fname=os.path.join(subdir, "%s.%s.svg"%(prefix, si)))
        else:
            pd_to_svg(pd, fname="%s.%s.svg"%(prefix, i))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Produce SVG images for all shadows in a PDSTOR file.")
    parser.add_argument('pdstor_f', type=argparse.FileType('r'),
                        help="a PDSTOR container to batch process")
    parser.add_argument('prefix', help="prefix for output SVG filenames")
    parser.add_argument('--skip_header', dest='read_header', action='store_const',
                        const=False, default=True, help="skip PDSTOR header")

    args = parser.parse_args()

    pdstor_to_svg(args.pdstor_f, args.read_header, prefix=args.prefix)
