from libpl.plcurve import *
import numpy as np
from numpy import linalg as LA
import pylab as P
import random

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

from scipy.stats import vonmises

class CrossingTetrahedron(object):
    def __init__(self, A, B, C, D):
        self._A = np.asarray(A)
        self._B = np.asarray(B)
        self._C = np.asarray(C)
        self._D = np.asarray(D)

    def _getset_dihedral(self, prop, A, B, C, D):
        ret = getattr(self, prop, None)
        if ret is None:
            ret = get_dihedral_angle(A,B,C,D)
            setattr(self, prop, ret)
        return ret

    def _getset_face(self, prop, A, B, C):
        ret = getattr(self, prop, None)
        if ret is None:
            ret = get_face_angle(A,B,C)
            setattr(self, prop, ret)
        return ret

    def _getset_edge(self, prop, A, B):
        ret = getattr(self, prop, None)
        if ret is None:
            ret = LA.norm( A-B )
            setattr(self, prop, ret)
        return ret

    def _getset_prop(self, prop, f, *args, **kwargs):
        def _setret():
            ret = f(*args, **kwargs)
            setattr(self, prop, ret)
            return ret
        return getattr(self, prop, _setret())

    @property
    def theta_0(self):
        return self._getset_dihedral("_theta_0", self._A, self._D, self._B, self._C)
    @property
    def theta_1(self):
        return self._getset_dihedral("_theta_1", self._A, self._B, self._C, self._D)
    @property
    def theta_2(self):
        return self._getset_dihedral("_theta_2", self._C, self._B, self._A, self._D)
    @property
    def theta_3(self):
        return self._getset_dihedral("_theta_3", self._C, self._D, self._A, self._B)

    @property
    def alpha_0(self):
        return self._getset_dihedral("_alpha_0", self._A, self._C, self._B, self._D)
    @property
    def alpha_1(self):
        return self._getset_dihedral("_alpha_1", self._B, self._D, self._A, self._C)

    @property
    def diao_0(self):
        return self._getset_face("_diao_0", self._D, self._A, self._C)
    @property
    def diao_1(self):
        return self._getset_face("_diao_1", self._B, self._D, self._A)

    @property
    def radius(self):
        return self._getset_edge("_radius", self._A, self._D)

    @classmethod
    def new_from_plc(cls, plc, i, k):
        verts = plc.components[0].vertices
        N = len(verts)

        A,C,B,D = (np.array(verts[i%N]), np.array(verts[(i+1)%N]),
                   np.array(verts[(i+k+1)%N]), np.array(verts[(i+k+2)%N]))

        del plc
        return cls(A,B,C,D)

def circmean(alpha, axis=None, ratio=1):
    mean_angle = np.arctan2(np.mean(np.sin(alpha*ratio),axis),
                            np.mean(np.cos(alpha*ratio),axis))
    return mean_angle/ratio

def get_dihedral_angle(A, B, C, D):
    """Return the dihedral angle from ABC to ABD"""
    v1 = C - A
    v1 = v1 / LA.norm(v1)
    v2 = B - A
    v2 = v2 / LA.norm(v2)
    v3 = D - A
    v3 = v3 / LA.norm(v3)

    n12 = np.cross(v1, v2)
    u12 = n12 / LA.norm(n12)

    n23 = np.cross(v2, v3)
    u23 = n23 / LA.norm(n23)

    return np.arccos(-np.dot(u12, u23))

def get_face_angle(A, B, C):
    """Returns angle <ABC."""
    u, v = A-B, C-B
    u, v = u/LA.norm(u), v/LA.norm(v)

    return np.arccos(-np.dot(u, v))

def random_tetrahedra_one_poly(n_edges, rng, num_results, i=None, k=None, min_r=None):
    rplc = PlCurve.random_equilateral_closed_polygon(n_edges, rng)
    for _ in range(num_results):
        _i = i if i is not None else random.randrange(n_edges)
        _k = k if k is not None else random.randrange(2,n_edges-3)
        yield CrossingTetrahedron.new_from_plc(rplc, _i, _k)

def random_tetrahedra(n_edges, rng, num_results, i=None, k=None, min_r=None):
    for _ in range(num_results):
        rplc = PlCurve.random_equilateral_closed_polygon(n_edges, rng)
        _i = i if i is not None else random.randrange(n_edges)
        _k = k if k is not None else random.randrange(2,n_edges-3)
        yield CrossingTetrahedron.new_from_plc(rplc, _i, _k)
        #del rplc

def hist3d(xvar, yvar, nboxes=10):
    hist, xedges, yedges = np.histogram2d(xvar, yvar, nboxes, normed=True)
    #print hist
    #hist = hist*1.0 / sum(hist)
    #print hist
    print hist.sum()/np.pi/np.pi

    elements = (len(xedges) - 1) * (len(yedges) - 1)
    xposm, yposm = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.25)

    xpos = xposm.flatten()
    ypos = yposm.flatten()
    zpos = np.zeros(elements)
    dx = 0.5*np.ones_like(zpos)
    dy = dx.copy()
    dz = hist.flatten()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #return ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
    #ax.scatter(xpos, ypos, dz)
    #ax.plot_trisurf(xpos, ypos, dz, cmap=cm.coolwarm)
    ax.plot_surface(xposm, yposm, hist, cstride=1, rstride=1)
    xhist = hist.sum(0)*nboxes/10
    xhist /= xhist.sum()/np.pi
    yhist = hist.sum(1)*nboxes/10
    yhist /= yhist.sum()/np.pi
    print xhist.sum()/np.pi
    print yhist.sum()/np.pi
    ax.plot(xpos[:nboxes], xhist, zdir="y", color="r")
    ax.plot(xpos[:nboxes], yhist, zdir="x", color="r")
    ax.set_xlim([0, np.pi])
    ax.set_ylim([0, np.pi])


if __name__ == "__main__":
    rng = RandomGenerator(444)
    theta = [[], [], [], []]
    diao_th = [[], []]
    EPS = 0.1
    ALPHA = np.pi/4
    for tet in random_tetrahedra(100, rng, 100000, min_r=4):
        #if t2 < ALPHA - EPS or t2 > ALPHA + EPS:
        #    continue
        theta[0].append(tet.theta_0)
        theta[1].append(tet.theta_1)
        #theta[2].append(tet.theta_2)
        theta[3].append(tet.theta_3)
        diao_th[0].append(tet.diao_0)
        diao_th[1].append(tet.diao_1)
        del tet
    theta[0] = np.array(theta[0])
    theta[1] = np.array(theta[1])
    theta[2] = np.array(theta[2])
    theta[3] = np.array(theta[3])

    #mirr_theta0 = np.concatenate((-theta[0], theta[0]))

    print "Done generating data"

    #PLOT = "VONMISES_THETA0"
    #PLOT = "THETA1 vs THETA3"
    PLOT = "FACE ANGLES"
    #PLOT = "PHI ANGLES"

    # Histogram of theta0 vs vonmises(mu=0, kappa=.5)

    if PLOT == "VONMISES_THETA0":
        n, bins, patches = P.hist(theta[0], 40, normed=1, histtype="stepfilled")
        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

        kappa = .5
        x = np.linspace(0,
                        vonmises.ppf(0.99, kappa), 100)
        rv = vonmises(kappa)

        P.plot(x, 2*rv.pdf(x), 'r-', lw=5, alpha=0.6, label="k = %d"%kappa)

        P.show()

    if PLOT == "FACE ANGLES":
        n, bins, patches = P.hist(diao_th[1], 20, normed=1, histtype="stepfilled")
        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        n, bins, patches = P.hist(diao_th[0], 20, normed=1, histtype="stepfilled")
        P.setp(patches, 'facecolor', 'b', 'alpha', 0.75)

        kappa = .5
        x = np.linspace(0,
                        vonmises.ppf(0.99, kappa), 100)
        rv = vonmises(kappa)

        P.plot(x, .5*np.sin(x), 'r-', lw=5, alpha=0.6, label="k = %d"%kappa)

        P.show()

    if PLOT == "PHI ANGLES":
        n, bins, patches = P.hist(theta[1], 40, normed=1, histtype="stepfilled")
        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        n, bins, patches = P.hist(theta[3], 40, normed=1, histtype="stepfilled")
        P.setp(patches, 'facecolor', 'b', 'alpha', 0.75)

        P.plot(1.0/np.pi, 'r-', lw=5, alpha=0.6)

        P.show()


    elif PLOT == "THETA1 vs THETA3":
        hist3d(theta[1], theta[3], nboxes=20)
        plt.show()

    #n, bins, patches = P.hist(theta[3], 20, normed=1, histtype="step")
    #P.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
    #n, bins, patches = P.hist(theta[0], 20, normed=1, histtype="step")
    #P.setp(patches, 'facecolor', 'm', 'alpha', 0.75)
    #n, bins, patches = P.hist(theta[2], 20, normed=1, histtype="step")
    #P.setp(patches, 'facecolor', 'm', 'alpha', 0.75)

    #y = P.plot(bins, [1.0/np.pi for _ in bins], 'r', linewidth=1.5)

    #plt.show()
