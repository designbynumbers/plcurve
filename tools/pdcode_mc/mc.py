import random

from libpl.pdcode import *

def shadow_r1_plus(pd, edge, face_pm):
    return pd.R1_loop_addition(*edge.face_index_pos()[face_pm])

def shadow_r1_minus(pd, edge, face_pm):
    if pd.ncross < 2:
        return pd

    face = edge.face_pos()[face_pm][0]
    if face.nedges == 1:
        loop_x = pd.crossings[face.vertices[0]]
        loop_x.sign = 0
        return pd.R1_loop_deletion(loop_x.index)

    return pd

def shadow_r2a_plus(pd, edge, face_pm):
    face, edge_1p = edge.face_pos()[face_pm]
    if face.nedges < 3:
        return pd

    face_delta = random.randint(1, len(face)-1)
    edge_2p = (edge_1p+face_delta)%len(face)

    print face, face.signs
    print edge_1p, edge_2p
    return pd.R2_bigon_addition(face.index, edge_1p, edge_2p, 2)

def shadow_r2a_minus(pd, edge, face_pm):
    if pd.ncross < 3:
        return pd
    face = edge.face_pos()[face_pm][0]

    print face, face.signs
    if len(face.vertices) == 2:
        return pd.R2_bigon_elimination_vertices(*face.vertices)[0]

    return pd

def shadow_r3(pd, edge, face_pm):
    return pd

    face = edge.face_pos()[face_pm][0]
    if face.nedges != 3:
        return pd

    tail_x = edge.prev_crossing()
    tail_p = edge.tailpos
    head_x = edge.next_crossing()
    head_p = edge.headpos

    opp_x, = [x for x in face.get_vertices() if not
              (x == tail_x or x == head_x)]

    outer_f, outer_p = edge.face_pos()[(face_pm+1)%2]

    if face_pm == PD_POS_ORIENTATION:
        e = edge.index
        b = head_x.edges[(head_p+1)%4]
        a = head_x.edges[(head_p+2)%4]
        f = head_x.edges[(head_p+3)%4]

        g = tail_x.edges[(tail_p+1)%4]
        d = tail_x.edges[(tail_p+2)%4]
        c = tail_x.edges[(tail_p+3)%4]


        print
    else:
        return pd

    pd.regenerate()
    return pd

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
        print transition

        self._diagram = transition(self._diagram, edge, face_pm)

state = PDMarkovState()
