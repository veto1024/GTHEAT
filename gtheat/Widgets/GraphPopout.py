#!/usr/bin/python

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import math
import numpy as np
from gtheat.Widgets.MplCanvas import MplCanvas

COLORS = ['r', 'g', 'blue', 'y', 'purple']
class GraphContainer(object):
    def __init__(self, x=None, y=None, title="", yAxisLabel="", xAxisLabel=r"$/rho$", legend=None):
        super(GraphContainer, self).__init__()
        if legend is None:
            self.legend = []
        else:
            self.legend = legend
        self.x = x
        self.y = y
        self.title = title
        self.yAxisLabel = yAxisLabel
        self.xAxisLabel = xAxisLabel

class InnerGraph(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.layout = QVBoxLayout()
        self.graph = MplCanvas()
        self.setup_Ui()

        if kwargs.get("title"):
            self.graph.set_window_title(kwargs.get("title"))


    def setup_Ui(self):
        self._create_spinboxes()
        self.layout.addWidget(self.graph)
        self.layout.addWidget(self.axesModifierFrame)
        self.setLayout(self.layout)
        self.setMinimumWidth(600)

    def _create_spinboxes(self):
        self.axesModifierLayout = QGridLayout()
        self.axesModifierFrame = QFrame()

        self.XAxisMin = QDoubleSpinBox()
        self.XAxisMin.setSingleStep(0.01)
        self.XAxisMin.setRange(0.0, 1.0)
        self.XAxisMin.setValue(0.85)

        self.XAxisMax = QDoubleSpinBox()
        self.XAxisMax.setRange(0.0, 1.0)
        self.XAxisMax.setSingleStep(0.01)
        self.XAxisMax.setValue(1.0)

        self.YAxisMax = QDoubleSpinBox()
        self.YAxisMax.setSingleStep(1.0)
        self.YAxisMax.setValue(10.0)
        self.YAxisMax.setMaximum(1E6)

        self.YAxisMin = QDoubleSpinBox()
        self.YAxisMin.setSingleStep(1.0)
        self.YAxisMin.setMinimum(-1E6)
        self.YAxisMin.setValue(0.0)

        self.XAxisLabel = QLabel("X-axis")
        self.YAxisLabel = QLabel("Y-axis")


        self.XAxisMin.valueChanged.connect(self.graph.updateXMin)
        self.XAxisMax.valueChanged.connect(self.graph.updateXMax)
        self.YAxisMin.valueChanged.connect(self.graph.updateYMin)
        self.YAxisMax.valueChanged.connect(self.graph.updateYMax)

        self.axesModifierLayout.addWidget(self.XAxisLabel, 0, 0)
        self.axesModifierLayout.addWidget(self.XAxisMin, 1, 0)
        self.axesModifierLayout.addWidget(self.XAxisMax, 1, 1)
        self.axesModifierLayout.addWidget(self.YAxisLabel, 2, 0)
        self.axesModifierLayout.addWidget(self.YAxisMin, 3, 0)
        self.axesModifierLayout.addWidget(self.YAxisMax, 3, 1)
        self.axesModifierLayout.setAlignment(Qt.AlignVCenter)
        self.axesModifierFrame.setLayout(self.axesModifierLayout)
        self.axesModifierFrame.setFrameStyle(QFrame.StyledPanel)
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(15)
        self.axesModifierFrame.setGraphicsEffect(shadow)
        self.axesModifierFrame.setMaximumWidth(200)

    def updateFig(self, d, **kwargs):
        maxY, oom = self._approxMax(np.max(d.y))
        self.graph.updateFig(d.x, d.y, legend=d.legend, **kwargs)
        self.graph.updateYMax(maxY)
        self.YAxisMax.setValue(maxY)
        self.YAxisMax.setSingleStep(10**(oom-1))
        self.YAxisMin.setSingleStep(10**(oom-1))
        self.graph.updateXMin(.85)
        print("Max approximated as:" + str(maxY))

    def _oom(self, number):
        return math.floor(math.log(number, 10))

    def _approxMax(self, number):
        oom = self._oom(number)
        if oom < 0:
            a = number // (10**oom)
        else:
            a = number // int(10**oom)
        return (a + 1) * int(10**oom), oom

    def set_title(self, title):
        self.graph.set_title(title)

class GraphPopout(QDialog):

    def __init__(self, data: [GraphContainer], grid=11, title="", parent=None):
        super(GraphPopout, self).__init__(parent)
        self.windowTitle = title
        self.gridSize = grid
        self.data = data
        self.setup_Ui()

    def setup_Ui(self):
        self.setWindowTitle(self.windowTitle)
        self.mainLayout = QVBoxLayout()
        self.graphLayout = QGridLayout()
        self._graphs = []


        if self.gridSize == 1:
            # self.shotDensityGraphs.updateFig(shot.rho, [shot.n.i.fsa, shot.n.e.fsa], legend=[r"$n_i$", r"$n_e$"], yFormatter=FormatStrFormatter('%.0E')
            self._graphs.append(InnerGraph())
            self.graphLayout.addWidget(self._graphs[0], 0, 0)
            self._graphs[0].updateFig(self.data[0])
        else:
            for n in range(self.gridSize):
                a = n // 2
                b = n % 2
                self._graphs.append(InnerGraph())
                self.graphLayout.addWidget(self._graphs[n], a, b)
                self._graphs[n].updateFig(self.data[n],
                                          keepLims=True,
                                          color=COLORS[n],)

        self.mainLayout.addLayout(self.graphLayout)
        self.setLayout(self.mainLayout)

