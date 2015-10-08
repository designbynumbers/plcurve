from __future__ import division
import random

from libpl.pdcode.diagram import *

MAX_DIA_N = 16

def shadow_r1_plus(pd, edge, face_pm, z, max_n=None):
    if max_n is not None:
        if pd.ncross + 1 > max_n:
            return pd, False

    alpha = random.random()
    if alpha > z:
        return pd, False

    return pd.R1_loop_addition(*edge.face_index_pos()[face_pm]), True

def shadow_r1_minus(pd, edge, face_pm, z, max_n=None):
    if pd.ncross < 2:
        return pd, False

    face = edge.face_pos()[face_pm][0]
    if face.nedges == 1:
        loop_x = pd.crossings[face.vertices[0]]
        loop_x.sign = 0
        return pd.R1_loop_deletion(loop_x.index), True

    return pd, False

def shadow_r2a_plus(pd, edge, face_pm, z, max_n=None):
    if max_n is not None:
        if pd.ncross+2 > max_n:
            return pd, False

    alpha = random.random()
    if alpha > z*z:
        return pd, False

    face, edge_1p = edge.face_pos()[face_pm]
    if face.nedges < 3:
        return pd, False

    face_delta = random.randint(1, len(face)-1)
    edge_2p = (edge_1p+face_delta)%len(face)

    #print face, face.signs
    #print edge_1p, edge_2p
    return pd.R2_bigon_addition(face.index, edge_1p, edge_2p, 2), True

def shadow_r2a_minus(pd, edge, face_pm, z, max_n=None):
    if pd.ncross < 3:
        return pd, False
    face = edge.face_pos()[face_pm][0]

    alpha = edge.prev_edge().face_pos()[(face_pm+1)%2][0].nedges
    beta = edge.next_edge().face_pos()[(face_pm+1)%2][0].nedges

    gamma = alpha + beta - 3
    if alpha == 1 or beta == 1 or gamma < 2:
        return pd, False
    kappa = random.random()
    if kappa > 1./gamma:
        return pd, False

    #print face, face.signs
    if len(face.vertices) == 2:
        return pd.R2_bigon_elimination_vertices(*face.vertices)[0], True

    return pd, False

def shadow_r3(pd, edge, face_pm, z, max_n=None):
    face = edge.face_pos()[face_pm][0]
    if face.nedges != 3 or len(set(face.vertices)) != 3:
        return pd, False

    #print repr(pd)
    #print edge, face, face.nedges
    return pd.R3_triangle_flip(face.index), True

class PDMarkovState(object):
    transitions = (
        shadow_r1_plus,
        shadow_r1_minus,
        shadow_r2a_plus,
        shadow_r2a_minus,
        shadow_r3,
    )

    def __init__(self, initial=None):
        if initial is None:
            initial = PlanarDiagram.unknot(1)
            initial.equiv_type = UO_ISOMORPHISM
        for x in initial.crossings:
            x.sign = 2
        self._diagram = initial

    def step(self, z=0.5, max_n=16):
        edge = random.choice(self._diagram.edges)
        face_pm = random.choice((0,1))
        #p = random.random()
        transition = random.choice(self.transitions)
        #print transition

        self._diagram, success = transition(self._diagram, edge, face_pm, z, max_n)
        self._diagram.equiv_type = UO_ISOMORPHISM
        return transition, success
        #print self._diagram.is_ok()

from collections import defaultdict

def run(z, N=200000, n_cross=2, mix_steps=100, max_n=16):
    state = PDMarkovState()
    MIX_N = mix_steps
    mix_n = 0

    hashes = defaultdict(lambda: 0)
    dts = dict()
    transitions = defaultdict(lambda: defaultdict(lambda: 0))

    for i in range(N):
        if (i+1)%10000 == 0:
            print i
        #assert(state._diagram.ncross <= MAX_DIA_N)
        #from_hash = state._diagram.hash
        state.step(z=z, max_n=max_n)
        state._diagram.regenerate_hash()
        #to_hash = state._diagram.hash
        #transitions[from_hash][to_hash] += 1
        #state._diagram.regenerate()
        #if (i+1)%mix_n == 0:
            #knot = state._diagram.copy()
            #knot.randomly_assign_crossings()
        if mix_n == MIX_N:
            #and state._diagram.ncross == n_cross:
            mix_n = 0

            if state._diagram.ncross == n_cross:
                print i+1, state._diagram.ncross
                H = state._diagram#.hash
                hashes[H] += 1
                if H not in dts:
                    dts[H] = state._diagram.copy()
                #initial = PlanarDiagram.unknot(1)
                #for x in initial.crossings:
                #    x.sign = 2
                #state._diagram = initial

        elif mix_n < MIX_N:
            mix_n += 1

        if state._diagram.ncross > 200:
            break

    print hashes
    for hsh, D in dts.iteritems():
        n_roots = n_cross*2*2/len(D.build_map_isomorphisms(D))
        print hsh, hashes[hsh], n_roots, hashes[hsh]/n_roots
        #D.edit_copy()
    return hashes, transitions
