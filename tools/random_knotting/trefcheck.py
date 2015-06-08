from libpl import pdcode
from collections import defaultdict
import signal
import sys

def end():
    for kt in sorted(knotflys.iterkeys()):
        print
        print kt
        print "n\tcount\ttotal\t\t%s"%kt
        for ncross in range(MIN_X, MAX_X+1):
            if counts[ncross] <= 0:
                continue
            print "%s\t%s\t%s\t\t%s" %(
                ncross, int(match[kt][ncross]), int(counts[ncross]),
                match[kt][ncross]/counts[ncross])

def terminate(signal, frame):
    end()
    sys.exit(0)
signal.signal(signal.SIGINT, terminate)

N_TRIES = 20000

A_0_1 = pdcode.HOMFLYPolynomial("1")
M_3_1 = pdcode.HOMFLYPolynomial("-a^{-4} + -2a^{-2} + a^{-2}z^{2}")
P_3_1 = pdcode.HOMFLYPolynomial("-2a^{2} + a^{2}z^{2} + -a^{4}")
A_4_1 = pdcode.HOMFLYPolynomial("-a^{-2} + -1 + z^{2} + -a^{2}")
P_5_1 = pdcode.HOMFLYPolynomial("3a^{4} + -4a^{4}z^{2} + a^{4}z^{4} + 2a^{6} + -a^{6}z^{2}")
M_5_1 = pdcode.HOMFLYPolynomial("2a^{-6} + -a^{-6}z^{2} + 3a^{-4} + -4a^{-4}z^{2} + a^{-4}z^{4}")
P_5_2 = pdcode.HOMFLYPolynomial("-a^{2} + a^{2}z^{2} + a^{4} + -a^{4}z^{2} + a^{6}")
M_5_2 = pdcode.HOMFLYPolynomial("a^{-6} + a^{-4} + -a^{-4}z^{2} + -a^{-2} + a^{-2}z^{2}")
P_3_1_3_1 = pdcode.HOMFLYPolynomial("4a^{4} + -4a^{4}z^{2} + a^{4}z^{4} + 4a^{6} + -2a^{6}z^{2} + a^{8}")
M_3_1_3_1 = pdcode.HOMFLYPolynomial("a^{-8} + 4a^{-6} + -2a^{-6}z^{2} + 4a^{-4} + -4a^{-4}z^{2} + a^{-4}z^{4}")
A_3_1_3_1M = pdcode.HOMFLYPolynomial("2a^{-2} + -a^{-2}z^{2} + 5 + -4z^{2} + z^{4} + 2a^{2} + -a^{2}z^{2}")

knotflys = {
    "0_1": (A_0_1,),
    "3_1": (P_3_1, M_3_1,),
    "4_1": (A_4_1,),
    "5_1": (P_5_1, M_5_1),
    "5_2": (P_5_2, M_5_2),
    "3_1#3_1": (P_3_1_3_1, M_3_1_3_1),
    "3_1#3_1*": (A_3_1_3_1M,)
}

counts = defaultdict(lambda:0.0)

match = defaultdict(lambda:defaultdict(lambda:0.0))

MIN_X = 3
MAX_X = 100
for ncross in range(MIN_X, MAX_X+1):
    print "Testing diagrams in %s"%ncross
    print counts[ncross]
    for i in range(N_TRIES):
        rnd_pd = pdcode.PlanarDiagram.random_diagram(ncross, n_components=1, max_att=150)
        if rnd_pd:
            counts[ncross] += 1
            hfly = rnd_pd.homfly()
            if hfly and not hfly == A_0_1:
                for ktname, homflys in knotflys.iteritems():
                    if hfly in homflys:
                        match[ktname][ncross] += 1
                        break
            elif hfly and hfly == A_0_1:
                match['0_1'][ncross] += 1

end()
