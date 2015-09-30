import random

from libpl.pdcode import *

def shadow_r1_plus(pd, edge, face_pm, z):
    alpha = random.random()
    if alpha > z:
        return pd

    return pd.R1_loop_addition(*edge.face_index_pos()[face_pm])

def shadow_r1_minus(pd, edge, face_pm, z):
    if pd.ncross < 2:
        return pd

    face = edge.face_pos()[face_pm][0]
    if face.nedges == 1:
        loop_x = pd.crossings[face.vertices[0]]
        loop_x.sign = 0
        return pd.R1_loop_deletion(loop_x.index)

    return pd

def shadow_r2a_plus(pd, edge, face_pm, z):
    alpha = random.random()
    if alpha > z*z:
        return pd

    face, edge_1p = edge.face_pos()[face_pm]
    if face.nedges < 2:
        return pd

    face_delta = random.randint(1, len(face)-1)
    edge_2p = (edge_1p+face_delta)%len(face)

    #print face, face.signs
    #print edge_1p, edge_2p
    return pd.R2_bigon_addition(face.index, edge_1p, edge_2p, 2)

def shadow_r2a_minus(pd, edge, face_pm, z):
    if pd.ncross < 3:
        return pd
    face = edge.face_pos()[face_pm][0]

    #print face, face.signs
    if len(face.vertices) == 2:
        return pd.R2_bigon_elimination_vertices(*face.vertices)[0]

    return pd

def shadow_r3(pd, edge, face_pm, z):
    face = edge.face_pos()[face_pm][0]
    if len(set(face.vertices)) != 3 or face.nedges != 3:
        return pd

    #print repr(pd)
    #print edge, face, face.nedges
    return pd.R3_triangle_flip(face.index)

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
        for x in initial.crossings:
            x.sign = 2
        self._diagram = initial

    def step(self, z=0.5):
        edge = random.choice(self._diagram.edges)
        face_pm = random.choice((0,1))
        p = random.random()
        transition = random.choice(self.transitions)
        #print transition

        self._diagram = transition(self._diagram, edge, face_pm, z)

state = PDMarkovState()
mix_n = 400
for i in range(40000):
    if not state._diagram.is_ok():
        break
    print state._diagram.ncross
    state.step(z=0.17)
    state._diagram.regenerate()
    knot = state._diagram.copy()
    knot.randomly_assign_crossings()
    print knot.homfly()
    #print repr(state._diagram)
