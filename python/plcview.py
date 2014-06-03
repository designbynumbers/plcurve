import sys
import math
from math import sin, cos

from PySide.QtCore import *
from PySide.QtGui import *
from PySide.QtOpenGL import *
from OpenGL.GL import *

from time import sleep
import numpy as np
from scipy.spatial import ConvexHull

import libplcurve.plcurve as pl

def normalize(v):
    norm=np.linalg.norm(v)
    if norm==0:
        return v
    return v/norm


class Polygon(object):
    def __init__(self, plc):
        self._plc = plc
        self.colors = [(1,0,0), (1,1,0), (0,1,0),
                       (0,1,1), (0,0,1), (1,0,1)]

    def draw(self):
        for component,color in zip(self._plc,self.colors):
            glBegin(GL_LINE_LOOP)
            glColor3f( *color )
            for (x,y,z) in component:
                glVertex3f( x,y,z )
            glEnd()

class GLHull(object):
    def __init__(self, hull):
        hull.points=hull.points*1.01
        self._hull = hull

    def draw(self):
        for face in self._hull.simplices:
            glBegin(GL_LINE_LOOP)
            glColor3f( *(1,1,1) )
            for point in self._hull.points[face]:
                glVertex3f( *point )
            glEnd()

    def drawNumbers(self, prnt):
        for vtx in self._hull.vertices:
            prnt(self._hull.points[vtx]*1.1, "%s"%vtx)

class Window(QWidget):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.glWidget = GLWidget()
        self.thetaSlider = self.createSlider(SIGNAL("thetaChanged(int)"),
                                             self.glWidget.setTheta)
        self.phiSlider = self.createSlider(SIGNAL("phiChanged(int)"),
                                           self.glWidget.setPhi)

        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.glWidget)
        mainLayout.addWidget(self.thetaSlider)
        mainLayout.addWidget(self.phiSlider)
        self.setLayout(mainLayout)

        self.thetaSlider.setValue(0)

        self.setWindowTitle(self.tr("Polygon visualizer"))

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

def rotation_mtx(axis, a):
    n = normalize(axis)
    return [[n[0]*n[0] + cos(a) * (1 - n[0]*n[0]),
            n[0]*n[1] * (1 - cos(a)) - n[2] * sin(a),
            n[0]*n[2] * (1 - cos(a)) + n[1] * sin(a),
            0],

            [n[0] * n[1] * (1 - cos(a)) + n[2] * sin(a),
            n[1] * n[1] + cos(a) * (1 - n[1] * n[1]),
            n[1] * n[2] * (1 - cos(a)) - n[0] * sin(a),
            0],

            [n[0] * n[2] * (1 - cos(a)) - n[1] * sin(a),
            n[1] * n[2] * (1 - cos(a)) + n[0] * sin(a),
            n[2] * n[2] + cos(a) * (1 - n[2]*n[2]),
            0],

            [0, 0, 0, 1]]

class GLWidget(QGLWidget):
    def __init__(self, parent=None):
        super(GLWidget, self).__init__(parent)
        self._oldtheta = 0
        self._oldphi = 0
        self.chord_num=1
        self._tracking = False
        self.side = 3
        self.update_cbs = []
        self.rmtx = None

    def initializeGL(self):
        glClearColor(0.0, 0.0, 0.0, 0.0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LINE_SMOOTH)
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    def resizeGL(self, width, height):
        side = min(width, height)
        glViewport((width - side) / 2, (height - side) / 2, side, side)
        self.side = side
        self.width = width
        self.height = height

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        # glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
        glMatrixMode(GL_MODELVIEW)
        pass

    def printText(self, pt, s):
        x,y,z = pt
        glDisable(GL_DEPTH_TEST)
        self.qglColor(Qt.white)
        self.renderText(x,y,z, s)
        glEnable(GL_DEPTH_TEST)


    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        for poly in self.polys:
            poly.draw()
        self.polys[1].drawNumbers(self.printText)

    def mousePressEvent(self, event):
        if event.buttons() & Qt.LeftButton:
            self._tracking = True
            self._track_start = event.pos()

    def _trackballVector(self, qpoint):
        _x = -(qpoint.x()*2.0 - self.width)/self.side
        _y = (qpoint.y()*2.0 - self.height)/self.side
        _rsq = _x**2 + _y**2

        COMPOSITE = True
        if COMPOSITE:
            if _rsq <= 1.0/2:
                _z = math.sqrt(1.0 - _rsq)
            else:
                _z = 1.0/(2.0*math.sqrt(_rsq))
        else:
            if _rsq > 1:
                _z = 0
            else:
                _z = math.sqrt(1.0 - _rsq)

        return normalize(np.array((_x,_y,_z)))

    def mouseMoveEvent(self, event):
        if (self._tracking and
            event.buttons() & Qt.LeftButton and
            event.pos() != self._track_start):

            # Get start point
            start = self._trackballVector(self._track_start)

            # Get end point
            stop = self._trackballVector(event.pos())

            if not (np.allclose(start, stop)):
                theta = math.acos(np.dot(start, stop))
                axis = np.cross(stop, start)
                rmtx = rotation_mtx(axis,
                                    theta/2.0)
                mtx = glGetFloatv(GL_MODELVIEW_MATRIX)
                glLoadIdentity()
                glMultMatrixf(rmtx)
                glMultMatrixf(mtx)
                self.updateGL()
                self._track_start = event.pos()

    def setTheta(self, theta):
        theta_delta = self._oldtheta-theta
        self._oldtheta = theta
        theta = theta_delta/100.0*2*math.pi
        glMultMatrixf([cos(theta), -sin(theta), 0, 0,
                      sin(theta), cos(theta), 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1])
        self.updateGL()

    def setPhi(self, phi):
        phi_delta = self._oldphi-phi
        self._oldphi = phi
        phi = phi_delta/100.0*2*math.pi
        glMultMatrixf([1, 0, 0, 0,
                      0, cos(phi), -sin(phi), 0,
                      0, sin(phi), cos(phi), 0,
                      0, 0, 0, 1])
        self.updateGL()

    def updateGL(self):
        super(GLWidget, self).updateGL()

        for cb in self.update_cbs:
            if cb:
                cb(glGetFloatv(GL_MODELVIEW_MATRIX))

class PlCurveViewer(object):
    def __init__(self, plc, app):
        self.app = app
        self.rho = Polygon(plc)
        self.ch = GLHull(ConvexHull(plc[0]))

    def run(self, cb):
        self.window = Window()
        self.window.resize(400,400)
        self.window.show()

        self.window.glWidget.update_cbs = [cb]
        self.window.glWidget.polys = [self.rho, self.ch]
        self.app.exec_()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = Window()
    window.resize(400,400)
    window.show()

    rng = pl.RandomGenerator()
    rng.set(401)
    plc = pl.PlCurve.random_closed_polygon(rng, 2**12)
    plc.scale(24)
    rho = Polygon(plc)
    window.glWidget.polys = [rho]

    app.exec_()

    sys.exit()
