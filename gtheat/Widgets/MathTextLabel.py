#!/usr/bin/python

from PyQt5.QtWidgets import QWidget, QVBoxLayout

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


class MathTextLabel(QWidget):

    def __init__(self, mathText, size, parent=None, **kwargs):
        QWidget.__init__(self, parent, **kwargs)

        l=QVBoxLayout(self)
        l.setContentsMargins(0,0,0,0)

        r,g,b,a=self.palette().base().color().getRgbF()

        self._figure=Figure(edgecolor=(r,g,b), facecolor=(r,g,b))
        self._canvas=FigureCanvas(self._figure)
        l.addWidget(self._canvas)

        self._figure.clear()
        text=self._figure.suptitle(
            mathText,
            x=0.0,
            y=1.0,
            horizontalalignment='left',
            verticalalignment='top',
            size=size)
        self._canvas.draw()

        (x0,y0),(x1,y1)=text.get_window_extent().get_points()
        w=x1-x0; h=y1-y0

        self._figure.set_size_inches(w/80, h/80)
        self.setFixedSize(w,h)