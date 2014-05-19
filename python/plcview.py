import sys
import math

from PySide.QtCore import *
from PySide.QtGui import *
from PySide.QtOpenGL import *
from OpenGL.GL import *

from time import sleep

import libplcurve.plcurve as pl

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

class Window(QWidget):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.glWidget = GLWidget()
        self.thetaSlider = self.createSlider(SIGNAL("thetaChanged(int)"),
                                             self.glWidget.setTheta)

        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.glWidget)
        mainLayout.addWidget(self.thetaSlider)
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

class GLWidget(QGLWidget):
    def __init__(self, parent=None):
        super(GLWidget, self).__init__(parent)
        self._oldtheta = 0
        self.chord_num=1

    def initializeGL(self):
        glClearColor(0.0, 0.0, 0.0, 0.0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LINE_SMOOTH)
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    def resizeGL(self, width, height):
        side = min(width, height)
        glViewport((width - side) / 2, (height - side) / 2, side, side)

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        # glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
        #glMatrixMode(GL_MODELVIEW)
        pass

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        for poly in self.polys:
            poly.draw()

    def mousePressEvent(self, event):
        for poly in self.polys:
            print poly._plc
            poly._plc.append(poly._plc[0][::2**self.chord_num])
            self.chord_num+=1
        self.updateGL()

    def setTheta(self, theta):
        theta_delta = self._oldtheta-theta
        self._oldtheta = theta
        for poly in self.polys:
            poly._plc.rotate((1/math.sqrt(2),1/math.sqrt(2),0),theta_delta*2*math.pi/100.0)
        self.updateGL()

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
