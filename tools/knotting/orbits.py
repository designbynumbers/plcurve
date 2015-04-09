from libpl.pdcode import *

def iter_flags(pd):
    for cross in pd.crossings:
        for edge in set(cross):
            yield cross.index, 0, edge.index, edge.face_index_pos()[0][0]
            yield cross.index, 0, edge.index, edge.face_index_pos()[1][0]
            yield cross.index, 1, edge.index, edge.face_index_pos()[0][0]
            yield cross.index, 1, edge.index, edge.face_index_pos()[1][0]

#def iter_flags(pd):
#    for cross in pd.crossings:
#        for edge in set(cross):
#            yield cross.index, edge.index, edge.face_index_pos()[0][0]
#            yield cross.index, edge.index, edge.face_index_pos()[1][0]

def all_flags(pd):
    return list(iter_flags(pd))

def flag_orbit(pd, flag):
    x, s, e, f = flag
    return set((
        g.cross_maps[x],
        (s+g.edge_maps_ori[e]+1)%2,
        g.edge_maps[e],
        g.face_maps[f]) for
               g in pd.build_automorphisms())

#def flag_orbit(pd, flag):
#    x, e, f = flag
#    return set((g.cross_maps[x], g.edge_maps[e], g.face_maps[f]) for
#               g in pd.build_automorphisms())


def iter_orbits(pd):
    flags = set(all_flags(pd))
    while flags:
        flag = flags.pop()
        orbit = flag_orbit(pd, flag)
        yield orbit
        flags -= orbit
