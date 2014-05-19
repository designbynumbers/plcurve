import sys
import math

from PySide.QtCore import *
from PySide.QtGui import *
from PySide.QtOpenGL import *
from OpenGL.GL import *

from time import sleep
import random

import libplcurve.plcurve as pl
import libplcurve.pd as pd

from math import sin,cos

import numpy as np
from numpy import linalg

from itertools import izip

class Polygon(object):
    def __init__(self, plc):
        self._plc = plc

    def draw(self):
        glBegin(GL_LINE_LOOP)
        for (x,y,z) in self._plc.components[0].vertices:
            glVertex3f( x,y,z )
        glEnd()

dirs = ((.25,0.0,0.0),
        (0.0,.25,0.0),
        (-.25,0.0,0.0),
        (0.0,-.25,0.0))

def rotation_mtx(theta):
    return np.array([[cos(theta), sin(theta), 0], [-sin(theta), cos(theta), 0], [0,0,1]])

def dist(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

class Diagram(object):
    def __init__(self, dia):
        self._pd = dia
        self.reloc_vtx()
        self.selected = None
        self._dirs = [dirs for vt in self.verts]
        self.draw_points = True

    def tail_over(self, e):
        v = e.tail
        p = e.tailpos
        np = (p+1)%4
        ne = self._pd.edges[self._pd.crossings[v][np]]
        return (ne.tail == v and ne.tailpos == np)

    def head_over(self, e):
        v = e.head
        p = e.headpos
        np = (p+1)%4
        ne = self._pd.edges[self._pd.crossings[v][np]]
        return (ne.head == v and ne.headpos == np)

    def drawArc(self, edge, rgb, N=20):
        s, sp, e, ep = (edge.tail, edge.tailpos, edge.head, edge.headpos)
        start, end = self.verts[s], self.verts[e]
        _s = (start[0],start[1],0.0)
        _e = (end[0],end[1],0.0)
        to, ho = self.tail_over(edge), self.head_over(edge)
        ctrl_pts = (_s,
                    tuple(a+b for a, b in zip(_s,self._dirs[s][sp])),
                    tuple(a+b for a, b in zip(_e,self._dirs[e][ep])),
                    _e)
        colors = []
        if to:
            colors.append(rgb + (0.0,))
        colors.append(rgb + (1.0,))
        if to or ho:
            for _ in range(5):
                colors.append(rgb + (1.0,))
        if ho:
            colors.append(rgb + (0.0,))
            print colors
        glMap1f(GL_MAP1_VERTEX_3,
                0.0, 1.0, ctrl_pts)
        glEnable(GL_MAP1_VERTEX_3)
        glMap1f(GL_MAP1_COLOR_4,
                0.0, 1.0,
                colors)
        glEnable(GL_MAP1_COLOR_4)
        glBegin(GL_LINE_STRIP)
        for i in range(N+1):
            glEvalCoord1f(i*1.0/N)
        glEnd()

    def reloc_vtx(self):
        _rndpt = lambda: (1.6*random.random()-.8,1.6*random.random()-.8)
        self.verts = [_rndpt() for _ in
                          range(dia.ncross)]

    def move_selected(self, dx, dy):
        if self.selected is not None:
            self.verts[self.selected] = (self.verts[self.selected][0]+dx,
                                         self.verts[self.selected][1]+dy)

    def set_rot(self, theta):
        if self.selected is None:
            return
        A = rotation_mtx(theta)
        self._dirs[self.selected] = tuple(np.dot(A, e) for e in dirs)

    def place_selected(self, x, y):
        if self.selected is not None:
            self.verts[self.selected] = (x,y)

    def select_near(self, x, y, delta):
        for i, vert in enumerate(self.verts):
            d = dist(vert, (x,y))
            if (d < delta):
                self.selected=i
                return
        self.selected=None

    def rademacher_I(self):
        if not self.selected: return
        cross = self._pd.crossings[self.selected]
        # Find loop which can be RI'd out
        guess = self._pd.edges[cross[0]]
        if guess.head != guess.tail:
            guess = self._pd.edges[cross[2]]
            if guess.head != guess.tail:
                return # no loop here

        # guess is loop to remove
        out_edge = (guess.headpos+2)%4
        in_edge = self._pd.edges[(guess.tailpos+2)%4]
        # we will join these two edges
        self._pd.del_crossing(self.selected)
        self._pd.crossings[in_edge.tail][in_edge.tailpos] = out_edge
        self._pd.regenerate()

    def draw(self):
        colors = ((1.0,0.0,0.0),
                  (0.0,1.0,0.0),
                  (0.0,0.0,1.0))
        i = 0
        for color, comp in izip(colors, self._pd.components):
            i+=1
            for edge in comp:
                edge = self._pd.edges[edge]
                self.drawArc(edge, color)
        if self.draw_points:
            glBegin(GL_POINTS)
            glColor3f(1.0,1.0,1.0)
            for i, vtx in enumerate(self.verts):
                if (self.selected == i):
                    glColor3f(1.0,0.0,0.0)
                    glVertex2f(*vtx)
                    glColor3f(1.0,1.0,1.0)
                else:
                    glVertex2f(*vtx)
            glEnd()


class Window(QWidget):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.glWidget = GLWidget()
        self.thetaSlider = self.createSlider(SIGNAL("thetaChanged(int)"),
                                             self.glWidget.setTheta)
        self.rademIButton = self.createButton(SIGNAL("rademacherI(int)"),
                                              self.glWidget.radI)

        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.glWidget)
        mainLayout.addWidget(self.thetaSlider)
        mainLayout.addWidget(self.rademIButton)
        self.setLayout(mainLayout)

        self.thetaSlider.setValue(0)

        self.setWindowTitle(self.tr("PlanarDigram visualizer"))

    def createSlider(self, changedSignal, setterSlot):
        slider = QSlider(Qt.Vertical)

        slider.setRange(0, 100)
        slider.setSingleStep(1)
        slider.setPageStep(10)
        slider.setTickInterval(10)
        slider.setTickPosition(QSlider.TicksRight)

        self.glWidget.connect(slider, SIGNAL("valueChanged(int)"),
                              setterSlot)
        self.connect(self.glWidget, changedSignal, slider, SLOT("setValue(int)"))

        return slider

    def createButton(self, changedSignal, setterSlot):
        slider = QPushButton("&Rademacher", self)

        self.glWidget.connect(slider, SIGNAL("valueChanged(int)"),
                              setterSlot)
        self.connect(self.glWidget, changedSignal, slider, SLOT("setValue(int)"))

        return slider

class GLWidget(QGLWidget):
    def __init__(self, parent=None):
        super(GLWidget, self).__init__(parent)
        self._oldtheta = 0
        self.lastPos = QPoint()
        self.width, self.height, self.side = 300,300,300

    def initializeGL(self):
        glClearColor(0.0, 0.0, 0.0, 0.0)
        glEnable(GL_DEPTH_TEST)
        #glEnable(GL_LINE_SMOOTH)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
        glPointSize(6.0)

    def pixel_to_unit(self, x, y):
        return (2.0*x/self.side-1.0*self.width/self.side,
                1.0*self.height/self.side-2.0*y/self.side)

    def mousePressEvent(self, event):
        self.dia.select_near(*self.pixel_to_unit(event.x(), event.y()),
                             delta=300*.05/self.side)
        self.lastPos = QPoint(event.pos())
        self.updateGL()

    def mouseMoveEvent(self, event):
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()

        if event.buttons() & Qt.LeftButton:
            #self.polys[0].move_selected(dx/300.0, dy/300.0)
            self.dia.place_selected(*self.pixel_to_unit(event.x(), event.y()))
        self.updateGL()

        self.lastPos = QPoint(event.pos())

    def resizeGL(self, width, height):
        side = min(width, height)
        self.width, self.height, self.side = width,height,side
        glViewport((width - side) / 2, (height - side) / 2, side, side)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        #glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
        #glMatrixMode(GL_MODELVIEW)
        pass

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.dia.draw()

    def setTheta(self, theta):
        self.dia.set_rot(theta*2*math.pi/100)
        self.updateGL()

    def radI(self):
        self.dia.rademacher_I()
        self.updateGL()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = Window()
    window.resize(400,400)
    window.show()

    # rng = pl.RandomGenerator()
    # rng.set(401)
    # plc = pl.PlCurve.random_closed_polygon(rng, 50)
    # plc.scale(2.5)
    # rho = Polygon(plc)
    #f = open("8.ex.pdcode")
    #dia = pd.PlanarDiagram.read(f)
    r = pl.RandomGenerator()
    plc = pl.PlCurve.random_closed_polygon(r, 25)
    dia = plc.as_pd(r)
    #f.close()
    rho = Diagram(dia)
    window.glWidget.dia = rho

    app.exec_()

    sys.exit()
