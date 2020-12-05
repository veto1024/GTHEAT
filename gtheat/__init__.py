#!/usr/bin/python


#!/usr/bin/python

# My python thing
import sys
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import os.path as path
import numpy as np
from gtheat.lib.chi_i import chi_i
from lib.chi_i import chi_i
import pyqtgraph as pg
from lib.Widgets.Color import Color
from waitingspinnerwidget import QtWaitingSpinner
from lib.Widgets.GTLoaderDialog import GT3loaderDialog
import GT3
import os
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=7, height=10, dpi=100, fig=None):
        """

        :type fig: Figure
        """
        if fig:
            self.fig = fig.axes.figure
        else:
            self.fig = Figure(figsize=(width, height), dpi=dpi)
            self.axes = self.fig.add_subplot(111)
            self.noDataText = self.fig.text(0.4, 0.5, "No Data", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1})
        super(MplCanvas, self).__init__(self.fig)
        self.colorList = ["red", "blue", "green", "black", "purple"]

    def updateFig(self, rho, yvals, legend=None, keepLims=None, xFormatter: FormatStrFormatter = None, yFormatter: FormatStrFormatter=None, color="black"):
        try:
            for txt in self.fig.texts:
                txt.set_visible(False)
        except:
            pass
        if keepLims:
            xlim = self.axes.get_xlim()
            ylim = self.axes.get_ylim()
        self.axes.cla()

        self.axes.set_xticklabels(self.axes.get_xticks(), size=8)
        self.axes.set_yticklabels(self.axes.get_yticks(), size=8)
        if xFormatter:
            self.axes.xaxis.set_major_formatter(xFormatter)
        else:
            self.axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        if yFormatter:
            self.axes.yaxis.set_major_formatter(yFormatter)
        else:
            self.axes.yaxis.set_major_formatter(FormatStrFormatter('%.2E'))
        if type(yvals) == list:
            for n, yval in enumerate(yvals):
                self.axes.scatter(rho, yval, s=10, color=self.colorList[n])
        else:
            self.axes.scatter(rho, yvals, s=10, color=color)

        if legend:
            self.fig.legend(legend, fontsize=8)

        # for c in fig.axes.get_children():
        #      if isinstance(c, matplotlib.collections.PathCollection):
        #          offsets = np.array(c.get_offsets()).T
        #          self.axes.scatter(offsets[0], offsets[1], s=3)
        #self.fig.tight_layout(pad=.9)
        if keepLims:
            self.axes.set_xlim(xlim)
            self.axes.set_ylim(ylim)
        self.fig.canvas.draw()
        self.draw_idle()

    def add_scatter(self, x, y, color):
        self.axes.scatter(x, y, color=color, s=10)
        self.fig.canvas.draw()
        self.draw_idle()

    def updateXMin(self, val):
        xmax = self.axes.get_xlim()[1]
        self.axes.set_xlim(val, xmax)
        self.draw_idle()

    def updateXMax(self, val):
        xmin = self.axes.get_xlim()[0]
        self.axes.set_xlim(xmin, val)
        self.draw_idle()

    def updateYMin(self, val):
        ymax = self.axes.get_ylim()[1]
        self.axes.set_ylim(val, ymax)
        self.draw_idle()

    def updateYMax(self, val):
        ymin = self.axes.get_ylim()[0]
        self.axes.set_ylim(ymin, val)
        self.draw_idle()

class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setWindowTitle("GTHEAT 0.0.1 - A Cool Program")
        self.setStyleSheet("background-color: white;")
        self.resize(1280, 768)
        self.setContentsMargins(3, 3, 3, 3)

        self._create_toolbar()

        self.mainGraph = MplCanvas(self, width=2, height=6, dpi=100)

        # Create all the layouts
        self.mainLayout = QVBoxLayout()
        self.mainGraphLayout = QFormLayout()
        self.smallGraphsLayout = QHBoxLayout()
        self.leftPanelLayout = QVBoxLayout()
        self.LeftPanelAndMainLayout = QHBoxLayout()
        self.axesModifierLayout = QGridLayout()
        self.axesModifierFrame = QFrame()
        self.checkboxFrame = QFrame()
        self.checkboxLayout = QGridLayout()
        self.chiExpFrame = QFrame()
        self.chiExpLayout = QVBoxLayout()
        self.leftPanelFrame = QFrame()
        self.checkboxAndAxesModifierLayout = QVBoxLayout()

        # Set main widget and layout
        self.mainGraphLayout.addWidget(self.mainGraph)
        self.mainWidget = QWidget()
        self.mainWidget.setLayout(self.mainLayout)
        self.setCentralWidget(self.mainWidget)


        # Set geometries
        self.mainGraphLayout.setGeometry(QRect(220, 159, 551, 501))
        self.checkboxFrame.setMaximumSize(250, 150)
        self.axesModifierFrame.setMaximumSize(250, 150)
        self.smallGraphsLayout.setGeometry(QRect(50, 30, 930, 20))

        # Set widgets
        self.setNoDataGraphs()
        self._create_checkboxes()
        self._create_axis_modifiers()
        self._create_exp_checkboxes()

        # Add layouts
        self.checkboxAndAxesModifierLayout.addWidget(self.axesModifierFrame)
        self.checkboxAndAxesModifierLayout.addWidget(self.chiExpFrame)
        self.checkboxAndAxesModifierLayout.addWidget(self.checkboxFrame)
        self.leftPanelFrame.setLayout(self.checkboxAndAxesModifierLayout)
        self.leftPanelLayout.addWidget(self.leftPanelFrame)
        self.LeftPanelAndMainLayout.addLayout(self.leftPanelLayout)
        self.LeftPanelAndMainLayout.addLayout(self.mainGraphLayout, stretch=1)
        self.mainLayout.addLayout(self.smallGraphsLayout)
        self.mainLayout.addLayout(self.LeftPanelAndMainLayout)

    def setNoDataGraphs(self):

        smallHeight = 1
        smallWidth = 3

        self.shotDensityGraphs = MplCanvas(self, width=smallWidth, height=smallHeight, dpi=100)
        self.smallGraphsLayout.addWidget(self.shotDensityGraphs)

        self.shotTempGraph = MplCanvas(self, width=smallWidth, height=smallHeight, dpi=100)
        self.smallGraphsLayout.addWidget(self.shotTempGraph)

        self.shotRadiusGraph = MplCanvas(self, width=smallWidth, height=smallHeight, dpi=100)
        self.smallGraphsLayout.addWidget(self.shotRadiusGraph)

        self.shotVelocityGraph = MplCanvas(self, width=smallWidth, height=smallHeight, dpi=100)
        self.smallGraphsLayout.addWidget(self.shotVelocityGraph)

        self.shotSomethingGraph = MplCanvas(self, width=smallWidth, height=smallHeight, dpi=100)
        self.smallGraphsLayout.addWidget(self.shotSomethingGraph)

    def setSmallGraphs(self, shot: chi_i):
        """

        :param shot: The chi_i shot data
        """

        plt.close()
        plt.ioff()


        self.shotDensityGraphs.updateFig(shot.rho, [shot.n.i, shot.n.e], legend=[r"$n_i$", r"$n_e$"], yFormatter=FormatStrFormatter('%.0E'))
        self.shotTempGraph.updateFig(shot.rho, [shot.T.i.ev, shot.T.e.ev], legend=[r"$T_i$", r"$T_e$"], yFormatter=FormatStrFormatter('%.0E'))
        self.shotRadiusGraph.updateFig(shot.rho, shot.gyro_rad_ion, yFormatter=FormatStrFormatter('%.0E'))
        self.shotVelocityGraph.updateFig(shot.rho, shot.vth, yFormatter=FormatStrFormatter('%.0E'))
        self.shotSomethingGraph.updateFig(shot.rho, shot.vth, yFormatter=FormatStrFormatter('%.0E'))

    def onOpenGT3Started(self, s):
        fileName = QFileDialog.getOpenFileName(self,
                                               str("Open File"), os.getcwd(), str("GT3 files (*)"))
        if path.exists(fileName[0]):
            self.gt3 = GT3loaderDialog(fileName, parent=self)  # type: QWidget

            self.chi_i = self.gt3.chi_i # type: chi_i
            self.shot = self.gt3.shot # type: gt3
            self.rho = self.gt3.rho
            self.setSmallGraphs(self.chi_i)
            self.mainGraph.updateFig(self.rho, self.shot.rtrans.chi.i.chi4, yFormatter=FormatStrFormatter('%.1f'))
            self.mainGraph.axes.set_xlim((0.85, 1.))
            self.setCentralWidget(self.mainWidget)
            self._enable_checkboxes(True)
            self._enable_axis_modifiers(True)
            self._enable_exp_checkboxes(True)
        else:
            pass

    def createSpinner(self):
        self.spinner = QtWaitingSpinner(self)
        self.spinner.setRoundness(70.0)
        self.spinner.setMinimumTrailOpacity(15.0)
        self.spinner.setTrailFadePercentage(70.0)
        self.spinner.setNumberOfLines(12)
        self.spinner.setLineLength(10)
        self.spinner.setLineWidth(5)
        self.spinner.setInnerRadius(10)
        self.spinner.setRevolutionsPerSecond(1)
        self.spinner.setColor(QColor(81, 4, 71))

    def createMainGraph(self):

        self.graphWidget = pg.PlotWidget()
        self.graphWidget.setBackground('w')
        self.graphWidget.setXRange(0.85, 1.0)

    def mainReplot(self, int):
        neo = self.neoCheckbox.isChecked()
        itg12 = self.itg12Checkbox.isChecked()
        itg32 = self.itg32Checkbox.isChecked()
        gyro = self.gyroBohmCheckbox.isChecked()
        DA = self.driftAlfenCheckbox.isChecked()
        sum = self.combinedCheckbox.isChecked()
        chi1 = self.chi1RadioButton.isChecked()
        chi2 = self.chi2RadioButton.isChecked()
        chi3 = self.chi3RadioButton.isChecked()
        chi4 = self.chi4RadioButton.isChecked()

        if chi1:
            self.mainGraph.updateFig(self.rho, self.shot.rtrans.chi.i.chi1, yFormatter=FormatStrFormatter('%.1f'), keepLims=True, color="black")
        if chi2:
            self.mainGraph.updateFig(self.rho, self.shot.rtrans.chi.i.chi2, yFormatter=FormatStrFormatter('%.1f'), keepLims=True, color="black")
        if chi3:
            self.mainGraph.updateFig(self.rho, self.shot.rtrans.chi.i.chi3, yFormatter=FormatStrFormatter('%.1f'), keepLims=True, color="black")
        if chi4:
            self.mainGraph.updateFig(self.rho, self.shot.rtrans.chi.i.chi4, yFormatter=FormatStrFormatter('%.1f'), keepLims=True, color="black")

        if sum:
            self.mainGraph.add_scatter(self.rho, self.chi_i.plot_chis_custom(neo=neo,
                                                                             bohm=gyro,
                                                                             itg12=itg12,
                                                                             itg32=itg32,
                                                                             DA=DA,
                                                                             sum=True,
                                                                             show=False), color="red")
        else:
            if neo:
                self.mainGraph.add_scatter(self.rho, self.chi_i.neo_chi, color="red")
            if itg12:
                self.mainGraph.add_scatter(self.rho, self.chi_i.chi_itg_TS_12, color="green")
            if itg32:
                self.mainGraph.add_scatter(self.rho, self.chi_i.chi_itg_TS_32, color="purple")
            if DA:
                self.mainGraph.add_scatter(self.rho, self.chi_i.chi_DA, color="yellow")
            if gyro:
                self.mainGraph.add_scatter(self.rho, self.chi_i.gyro_bohm, color="blue")

    def _create_axis_modifiers(self):

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

        self.YAxisMin = QDoubleSpinBox()
        self.YAxisMin.setRange(-100., 100.)
        self.YAxisMin.setSingleStep(1.0)
        self.YAxisMin.setValue(0.0)

        self.XAxisLabel = QLabel("X-axis")
        self.YAxisLabel = QLabel("Y-axis")

        self.XAxisMin.valueChanged.connect(self.mainGraph.updateXMin)
        self.XAxisMax.valueChanged.connect(self.mainGraph.updateXMax)
        self.YAxisMin.valueChanged.connect(self.mainGraph.updateYMin)
        self.YAxisMax.valueChanged.connect(self.mainGraph.updateYMax)

        self.axesModifierLayout.addWidget(self.XAxisLabel, 0, 0)
        self.axesModifierLayout.addWidget(self.XAxisMin, 1, 0)
        self.axesModifierLayout.addWidget(self.XAxisMax, 1, 1)
        self.axesModifierLayout.addWidget(self.YAxisLabel, 2, 0)
        self.axesModifierLayout.addWidget(self.YAxisMin, 3, 0)
        self.axesModifierLayout.addWidget(self.YAxisMax, 3, 1)
        self.axesModifierLayout.setAlignment(Qt.AlignVCenter)
        self.axesModifierFrame.setLayout(self.axesModifierLayout)
        self.axesModifierFrame.setFrameStyle(QFrame.Panel)

        self._enable_axis_modifiers(False)

    def _enable_checkboxes(self, en=True):
        self.neoCheckbox.setEnabled(en)
        self.gyroBohmCheckbox.setEnabled(en)
        self.itg12Checkbox.setEnabled(en)
        self.itg32Checkbox.setEnabled(en)
        self.driftAlfenCheckbox.setEnabled(en)
        self.combinedCheckbox.setEnabled(en)

    def _enable_exp_checkboxes(self, en=True):
        self.chi1RadioButton.setEnabled(en)
        self.chi2RadioButton.setEnabled(en)
        self.chi3RadioButton.setEnabled(en)
        self.chi4RadioButton.setEnabled(en)

    def _enable_axis_modifiers(self, en=True):
        self.XAxisMin.setEnabled(en)
        self.XAxisMax.setEnabled(en)
        self.YAxisMin.setEnabled(en)
        self.YAxisMax.setEnabled(en)

    def _create_checkboxes(self):

        self.combinedCheckbox = QCheckBox(text="Sum Exp. Chis")
        self.neoCheckbox = QCheckBox(text="Neoclassical")
        self.neoCheckbox.setStyleSheet("color: red")

        self.gyroBohmCheckbox = QCheckBox(text="Gyro-Bohm")
        self.gyroBohmCheckbox.setStyleSheet("color: blue")

        self.itg12Checkbox = QCheckBox(text="ITG (1/2)")
        self.itg12Checkbox.setStyleSheet("color: green")

        self.itg32Checkbox = QCheckBox(text="ITG (3/2)")
        self.itg32Checkbox.setStyleSheet("color:purple")

        self.driftAlfenCheckbox = QCheckBox(text="Drift-Alfven")
        self.driftAlfenCheckbox.setStyleSheet("color: yellow")

        self._enable_checkboxes(False)

        self.combinedCheckbox.stateChanged.connect(self.mainReplot)
        self.neoCheckbox.stateChanged.connect(self.mainReplot)
        self.itg32Checkbox.stateChanged.connect(self.mainReplot)
        self.itg12Checkbox.stateChanged.connect(self.mainReplot)
        self.gyroBohmCheckbox.stateChanged.connect(self.mainReplot)
        self.driftAlfenCheckbox.stateChanged.connect(self.mainReplot)

        self.checkboxLayout.addWidget(self.neoCheckbox, 0, 0)
        self.checkboxLayout.addWidget(self.gyroBohmCheckbox, 0, 1)
        self.checkboxLayout.addWidget(self.itg12Checkbox, 1, 0)
        self.checkboxLayout.addWidget(self.itg32Checkbox, 1, 1)
        self.checkboxLayout.addWidget(self.driftAlfenCheckbox, 2, 0)
        self.checkboxLayout.addWidget(self.combinedCheckbox, 2, 1)
        self.checkboxLayout.setAlignment(Qt.AlignVCenter)


        self.checkboxFrame.setLayout(self.checkboxLayout)
        self.checkboxFrame.setFrameStyle(QFrame.Panel)

    def _create_exp_checkboxes(self):
        self.chi1RadioButton = QRadioButton(r"$q^{cond} = Q^{tot}$")
        self.chi2RadioButton = QRadioButton(r"$q^{cond} = Q^{tot} - Q^{conv}$")
        self.chi3RadioButton = QRadioButton(r"$q^{cond} = Q^{tot} - Q^{conv} - Q^{heatin}$")
        self.chi4RadioButton = QRadioButton(r"$q^{cond} = Q^{tot} - Q^{conv} - Q^{heatin} - Q^{visc}$")
        self.chi4RadioButton.setChecked(True)
        self._enable_exp_checkboxes(False)

        self.chi1RadioButton.toggled.connect(self.mainReplot)
        self.chi2RadioButton.toggled.connect(self.mainReplot)
        self.chi3RadioButton.toggled.connect(self.mainReplot)
        self.chi4RadioButton.toggled.connect(self.mainReplot)


        self.chiExpLayout.addWidget(self.chi1RadioButton)
        self.chiExpLayout.addWidget(self.chi2RadioButton)
        self.chiExpLayout.addWidget(self.chi3RadioButton)
        self.chiExpLayout.addWidget(self.chi4RadioButton)
        self.chiExpLayout.setAlignment(Qt.AlignVCenter)


        self.chiExpFrame.setLayout(self.chiExpLayout)
        self.chiExpFrame.setFrameStyle(QFrame.Panel)

    def _create_toolbar(self):
        # Set the toolbar
        self.toolbar = QToolBar("The Main Toolbar")

        #Set icon size
        self.toolbar.setIconSize(QSize(16, 16))

        # Create the toolbar
        self.toolbar = QToolBar("GTHEAT Toolbar")

        # Set icon size
        self.toolbar.setIconSize(QSize(16, 16))
        self.addToolBar(self.toolbar)

        icon = QIcon("gtheat/lib/icons/bug.png")
        open_GT3_action = QAction("Open GT3 File", self)
        open_GT3_action.setIcon(icon)
        open_GT3_action.setStatusTip("Open a GT3 File")
        open_GT3_action.triggered.connect(self.onOpenGT3Started)

        exitAct = QAction(QIcon('exit.png'), '&Exit', self)
        exitAct.setShortcut('Ctrl+Q')
        exitAct.setStatusTip('Exit application')
        exitAct.triggered.connect(qApp.quit)

        self.toolbar.addAction(open_GT3_action)

        menu = self.menuBar()
        menu.setNativeMenuBar(False)

        file_menu = menu.addMenu(u"&File")
        file_menu.addAction(open_GT3_action)
        file_menu.addAction(exitAct)


app = QApplication(sys.argv)
window = MainWindow()
window.show()
sys.exit(app.exec_())
