from libplcurve.pd import PlanarDiagram
from itertools import combinations, chain, compress, product
import pprint
pp = pprint.PrettyPrinter(indent=2)

def crossing_combinations(pdc):
    # Generate all possible crossing sign permutations.
    # First crossing is fixed as (+) so as to avoid mirror images
    for cross_set in product([0,1], repeat=len(pdc.crossings[1:])):
        new_pd = pdc.copy()
        for crs in compress(new_pd.crossings[1:], cross_set):
            crs.sign = (crs.sign+1)%2
        yield new_pd, cross_set

# TODO: All possible permutations
def component_combinations(pdc):
    # Generate all possible component permutations.
    for comp_set in product([0,1], repeat=len(pdc.components[1:])):
        new_pd = pdc.copy()
        for comp in compress(range(len(new_pd.components[1:])), comp_set):
            new_pd.reorient_component(comp, 0)
        yield new_pd, comp_set

def bin_list_to_int(blist):
    return sum((i*2)**n for i,n in enumerate(reversed(blist)))

def uid(pd, pos_in_stor, crs_sgn, cmp_sgn):
    return "_".join([str(len(pd.crossings)),
                     str(pos_in_stor),
                     str(bin_list_to_int(crs_sgn)),
                     str(len(pd.components)),
                     str(bin_list_to_int(cmp_sgn))])

def read_pdstor(f):
    # TODO: Actually parse the header
    header = f.readline()
    if "pdstor" not in header:
        return None # something bad actually happened
    # try to parse the data line
    header = f.readline().split()
    try:
        claimed, actual = header[1].split("/")
        nhashes = header[4]
    except:
        return None # actually throw an error here too
    print(claimed,actual,nhashes)

    # Read in the pd codes one-by-one.
    pd_dict = {}
    pd_id = 0
    pd = PlanarDiagram.read(f)
    pd.regenerate_hash()
    pd.printf("%PD")
    while pd is not None:
        for crs_pd, crs_sgn in crossing_combinations(pd):
            for cmp_pd, cmp_sgn in component_combinations(crs_pd):
                hsh = cmp_pd.homfly()
                print repr(cmp_pd)
                if hsh in pd_dict:
                    pd_dict[hsh].append(uid(cmp_pd, pd_id, crs_sgn, cmp_sgn))
                else:
                    pd_dict[hsh] = [uid(cmp_pd, pd_id, crs_sgn, cmp_sgn)]
        pd_id += 1
        pd = PlanarDiagram.read(f)
    return pd_dict

if __name__ == "__main__":
    f = file("../data/pdstors/3.pdstor")
    pds = read_pdstor(f)
    f.close()

    pp.pprint(pds)
