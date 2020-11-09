#!/usr/bin/python

from math import ceil

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

class Spinner(QWidget):
    def __init__(self, parent=None):
        super(Spinner, self).__init__(parent)
        #self.setAlignment(Qt.AlignCenter)
        self.pixmap = QPixmap("../icons/bug.png")

        self.setFixedSize(30, 30)
        self._angle = 0

        self.animation = QPropertyAnimation(self, b"angle", self)
        self.animation.setStartValue(0)
        self.animation.setEndValue(360)
        self.animation.setLoopCount(-1)
        self.animation.setDuration(2000)
        self.animation.start()


    @pyqtProperty(int)
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, value):
        self._angle = value
        self.update()


    def paintEvent(self, ev=None):
        painter = QPainter(self)
        painter.translate(15, 15)
        painter.rotate(self._angle)
        painter.translate(-15, -15)
        painter.drawPixmap(5, 5, self.pixmap)