#!/usr/bin/python

import sys
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import os.path as path
import numpy as np
from gtheat.lib.chi_i import chi_i
import pyqtgraph as pg
from gtheat.Widgets.GraphPopout import GraphContainer
from gtheat.Widgets.MplCanvas import MplCanvas
from gtheat.Widgets.Dialog import PreferencesDialog
from waitingspinnerwidget import QtWaitingSpinner
from gtheat.Widgets.GTLoaderDialog import GT3loaderDialog
from gtheat.Widgets.MathTextLabel import MathTextLabel
from gtheat.Widgets.GraphPopout import GraphPopout
import os
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')
from matplotlib.ticker import FormatStrFormatter


class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setWindowTitle("GTHEAT 0.2 - A Cool Program")
        self.setStyleSheet("background-color: white;")
        self.resize(1280, 768)
        self.setContentsMargins(3, 3, 3, 3)
        self.setStatusBar(QStatusBar())
        self.statusBar().setStyleSheet("border-top: 1px solid black")
        self.settings = QSettings("GT FRC", "GTHEAT")
        self._markerSize = 14
        self.mainGraph = MplCanvas(self, width=2, height=6, dpi=100)

        self._create_toolbar()

        # Create all the layouts
        self.mainLayout = QVBoxLayout()
        self.mainGraphLayout = QFormLayout()
        self.smallGraphsLayout = QHBoxLayout()
        self.leftPanelLayout = QVBoxLayout()
        self.LeftPanelAndMainLayout = QHBoxLayout()
        self.axesModifierLayout = QGridLayout()
        self.heatVisLayout = QGridLayout()
        self.heatVisFrame = QFrame()
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
        self._setNoDataGraphs()
        self._create_visc_modifiers()
        self._create_checkboxes()
        self._create_axis_modifiers()
        self._create_exp_checkboxes()

        # Add layouts
        self.checkboxAndAxesModifierLayout.addWidget(self.axesModifierFrame)
        self.checkboxAndAxesModifierLayout.addWidget(self.heatVisFrame)
        self.checkboxAndAxesModifierLayout.addWidget(self.chiExpFrame)
        self.checkboxAndAxesModifierLayout.addWidget(self.checkboxFrame)
        self.leftPanelFrame.setLayout(self.checkboxAndAxesModifierLayout)
        self.leftPanelLayout.addWidget(self.leftPanelFrame)
        self.LeftPanelAndMainLayout.addLayout(self.leftPanelLayout)
        self.LeftPanelAndMainLayout.addLayout(self.mainGraphLayout, stretch=1)
        self.mainLayout.addLayout(self.smallGraphsLayout)
        self.mainLayout.addLayout(self.LeftPanelAndMainLayout)

    def _setNoDataGraphs(self):

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

    def _setSmallGraphs(self, shot: chi_i):
        """

        :param shot: The chi_i shot data
        """

        plt.close()
        plt.ioff()


        self.shotDensityGraphs.updateFig(shot.rho, [shot.n.i.fsa, shot.n.e.fsa], yFormatter=FormatStrFormatter('%.0E'))
        self.shotTempGraph.updateFig(shot.rho, [shot.T.i.ev.fsa, shot.T.e.ev.fsa], yFormatter=FormatStrFormatter('%.0E'))
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
            self._setSmallGraphs(self.chi_i)
            self.mainGraph.updateFig(self.rho, self.shot.rtrans.chi.i.chi4, yFormatter=FormatStrFormatter('%.1f'))
            xAxisMin = self.XAxisMin.value()
            xAxisMax = self.XAxisMax.value()
            yAxisMin = self.YAxisMin.value()
            yAxisMax = self.YAxisMax.value()
            self.mainGraph.axes.set_xlim((xAxisMin, xAxisMax))
            self.mainGraph.axes.set_ylim((yAxisMin, yAxisMax))
            self.setCentralWidget(self.mainWidget)
            self._enable_checkboxes(True)
            self._enable_axis_modifiers(True)
            self._enable_exp_checkboxes(True)
            self.save_fig_action.setEnabled(True)
            # Enable Shot Data file menu
            self.open_velocity_graphs_action.setEnabled(True)
            self.open_temperature_graphs_action.setEnabled(True)
            self.open_scale_length_graphs_action.setEnabled(True)
            self.open_flux_graphs_action.setEnabled(True)
            self.open_electric_field_action.setEnabled(True)
            self.open_collisionality_action.setEnabled(True)
            self.open_density_graph_action.setEnabled(True)
        else:
            pass

    def _createSpinner(self):
        self.spinner = QtWaitingSpinner(parent=self)
        self.spinner.setRoundness(70.0)
        self.spinner.setMinimumTrailOpacity(15.0)
        self.spinner.setTrailFadePercentage(70.0)
        self.spinner.setNumberOfLines(12)
        self.spinner.setLineLength(10)
        self.spinner.setLineWidth(5)
        self.spinner.setInnerRadius(10)
        self.spinner.setRevolutionsPerSecond(1)
        self.spinner.setColor(QColor(81, 4, 71))

    # def createMainGraph(self):
    #
    #     self.graphWidget = pg.PlotWidget()
    #     self.graphWidget.setBackground('w')
    #     self.graphWidget.setXRange(0.85, 1.0)

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
        chiAll = self.chiShowAllButton.isChecked()
        viscTor = self.heatVisTor.value()
        viscPol = self.heatVisPol.value()
        viscChi = self.shot.rtrans._calc_chi_i_visc(vpolS=viscPol, vtorS=viscTor)
        s = self._markerSize
        if chiAll:
            self.mainGraph.updateFig(self.rho, None, color="black", keepLims=True, yFormatter=FormatStrFormatter('%.1f'))
            self.mainGraph.add_scatter(self.rho, self.shot.rtrans.chi.i.chi1, color="red", marker="x", s=s)
            self.mainGraph.add_scatter(self.rho, self.shot.rtrans.chi.i.chi2, color="green", marker="X", s=s)
            self.mainGraph.add_scatter(self.rho, self.shot.rtrans.chi.i.chi3, color="blue", marker="^", s=s)
            self.mainGraph.add_scatter(self.rho, viscChi, color="black", marker="2", s=s)
        else:
            if chi1:
                self.mainGraph.updateFig(self.rho, self.shot.rtrans.chi.i.chi1, yFormatter=FormatStrFormatter('%.1f'),
                                         keepLims=True, color="red", marker="x", s=s)
            if chi2:
                self.mainGraph.updateFig(self.rho, self.shot.rtrans.chi.i.chi2, yFormatter=FormatStrFormatter('%.1f'),
                                         keepLims=True, color="green", marker="X", s=s)
            if chi3:
                self.mainGraph.updateFig(self.rho, self.shot.rtrans.chi.i.chi3, yFormatter=FormatStrFormatter('%.1f'),
                                         keepLims=True, color="blue", marker="^", s=s)
            if chi4:
                self.mainGraph.updateFig(self.rho, viscChi, yFormatter=FormatStrFormatter('%.1f'),
                                         keepLims=True, color="black", marker="2", s=s)

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
                self.mainGraph.add_scatter(self.rho, self.chi_i.neo_chi, color="red", s=s)
            if itg12:
                self.mainGraph.add_scatter(self.rho, self.chi_i.chi_itg_TS_12, color="green", s=s)
            if itg32:
                self.mainGraph.add_scatter(self.rho, self.chi_i.chi_itg_TS_32, color="purple", s=s)
            if DA:
                self.mainGraph.add_scatter(self.rho, self.chi_i.chi_DA, color="yellow", s=s)
            if gyro:
                self.mainGraph.add_scatter(self.rho, self.chi_i.gyro_bohm, color="blue", s=s)

    def _create_visc_modifiers(self):

        self.heatVisTor = QDoubleSpinBox()
        self.heatVisPol = QDoubleSpinBox()

        self.heatVisTorLabel = MathTextLabel(r"$V_{\phi}^s$", 14)
        self.heatVisPolLabel = MathTextLabel(r"$V_{\theta}^s$", 14)

        self.heatVisTor.setSingleStep(0.01)
        self.heatVisTor.setRange(-1.0, 1.0)
        self.heatVisTor.setValue(0.1)

        self.heatVisPol.setSingleStep(0.01)
        self.heatVisPol.setRange(-1.0, 1.0)
        self.heatVisPol.setValue(0.1)

        self.heatVisTor.valueChanged.connect(self.mainReplot)
        self.heatVisPol.valueChanged.connect(self.mainReplot)


        self.heatVisLayout.addWidget(self.heatVisTorLabel, 0, 0)
        self.heatVisLayout.addWidget(self.heatVisPolLabel, 0, 1)
        self.heatVisLayout.addWidget(self.heatVisTor, 1, 0)
        self.heatVisLayout.addWidget(self.heatVisPol, 1, 1)
        self.heatVisFrame.setLayout(self.heatVisLayout)
        #self.heatVisLayout.setAlignment(Qt.AlignVCenter)
        self.heatVisFrame.setFrameStyle(QFrame.StyledPanel)
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(15)
        self.heatVisFrame.setGraphicsEffect(shadow)

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
        self.YAxisMax.setMaximum(1000.)

        self.YAxisMin = QDoubleSpinBox()
        self.YAxisMin.setSingleStep(1.0)
        self.YAxisMin.setMinimum(-100.)
        self.YAxisMin.setValue(0.0)

        self.LogScale = QCheckBox()
        self.LogScale.setChecked(False)

        self.XAxisLabel = QLabel("X-axis")
        self.YAxisLabel = QLabel("Y-axis")
        self.LogAxisLabel = QLabel("Log Scale")

        self.XAxisMin.valueChanged.connect(self.mainGraph.updateXMin)
        self.XAxisMax.valueChanged.connect(self.mainGraph.updateXMax)
        self.YAxisMin.valueChanged.connect(self.mainGraph.updateYMin)
        self.YAxisMax.valueChanged.connect(self.mainGraph.updateYMax)
        self.LogScale.stateChanged.connect(self.mainGraph.set_logPlot)

        self.axesModifierLayout.addWidget(self.XAxisLabel, 0, 0)
        self.axesModifierLayout.addWidget(self.XAxisMin, 1, 0)
        self.axesModifierLayout.addWidget(self.XAxisMax, 1, 1)
        self.axesModifierLayout.addWidget(self.YAxisLabel, 2, 0)
        self.axesModifierLayout.addWidget(self.YAxisMin, 3, 0)
        self.axesModifierLayout.addWidget(self.YAxisMax, 3, 1)
        self.axesModifierLayout.addWidget(self.LogAxisLabel, 0, 2)
        self.axesModifierLayout.addWidget(self.LogScale, 1, 2)
        self.axesModifierLayout.setAlignment(Qt.AlignVCenter)
        self.axesModifierFrame.setLayout(self.axesModifierLayout)
        self.axesModifierFrame.setFrameStyle(QFrame.StyledPanel)
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(15)
        self.axesModifierFrame.setGraphicsEffect(shadow)

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
        self.chiShowAllButton.setEnabled(en)

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
        self.checkboxFrame.setFrameStyle(QFrame.StyledPanel)
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(15)
        self.checkboxFrame.setGraphicsEffect(shadow)

    def _create_exp_checkboxes(self):
        self.chi1RadioButton = QRadioButton(r"$q^{cond} = Q^{tot}$")
        self.chi2RadioButton = QRadioButton(r"$q^{cond} = Q^{tot} - Q^{conv}$")
        self.chi3RadioButton = QRadioButton(r"$q^{cond} = Q^{tot} - Q^{conv} - Q^{heatin}$")
        self.chi4RadioButton = QRadioButton(r"$q^{cond} = Q^{tot} - Q^{conv} - Q^{heatin} - Q^{visc}$")
        self.chiShowAllButton = QRadioButton("Show All")
        self.chiShowAllButton.setChecked(True)

        # self.chi1RadioButton = QCheckBox(r"$q^{cond} = Q^{tot}$")
        # self.chi2RadioButton = QCheckBox(r"$q^{cond} = Q^{tot} - Q^{conv}$")
        # self.chi3RadioButton = QCheckBox(r"$q^{cond} = Q^{tot} - Q^{conv} - Q^{heatin}$")
        # self.chi4RadioButton = QCheckBox(r"$q^{cond} = Q^{tot} - Q^{conv} - Q^{heatin} - Q^{visc}$")
        # self.chiShowAllButton = QCheckBox("Show All")

        self.chi1RadioButton.setStyleSheet("color: red")
        self.chi2RadioButton.setStyleSheet("color: green")
        self.chi3RadioButton.setStyleSheet("color: blue")
        self.chi4RadioButton.setStyleSheet("color: black")
        self._enable_exp_checkboxes(False)


        self.chi1RadioButton.toggled.connect(self.mainReplot)
        self.chi2RadioButton.toggled.connect(self.mainReplot)
        self.chi3RadioButton.toggled.connect(self.mainReplot)
        self.chi4RadioButton.toggled.connect(self.mainReplot)
        self.chiShowAllButton.toggled.connect(self.mainReplot)


        self.chiExpLayout.addWidget(self.chi1RadioButton)
        self.chiExpLayout.addWidget(self.chi2RadioButton)
        self.chiExpLayout.addWidget(self.chi3RadioButton)
        self.chiExpLayout.addWidget(self.chi4RadioButton)
        self.chiExpLayout.addWidget(self.chiShowAllButton)
        self.chiExpLayout.setAlignment(Qt.AlignVCenter)


        self.chiExpFrame.setLayout(self.chiExpLayout)
        self.chiExpFrame.setFrameStyle(QFrame.StyledPanel)
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(15)
        self.chiExpFrame.setGraphicsEffect(shadow)

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

        # Create actions
        self.open_GT3_action = QAction("Open GT3 File", self)
        self.open_GT3_action.setIcon(self.style().standardIcon(QStyle.SP_DirOpenIcon))
        self.open_GT3_action.setStatusTip("Open a GT3 File")
        self.open_GT3_action.triggered.connect(self.onOpenGT3Started)

        self.save_fig_action = QAction("Save Figure", self)
        self.save_fig_action.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_fig_action.setStatusTip("Save Figure")
        self.save_fig_action.triggered.connect(self.saveFigureDialog)
        self.save_fig_action.setEnabled(False)

        self.set_preferences_action = QAction("Set Preferences", self)
        self.set_preferences_action.setStatusTip("Set Preferences")
        self.set_preferences_action.triggered.connect(self.set_preferences)

        self.open_velocity_graphs_action = QAction("Velocity Graphs", self)
        self.open_velocity_graphs_action.setStatusTip("Velocity Graphs")
        self.open_velocity_graphs_action.triggered.connect(self.open_velocity_graphs)
        self.open_velocity_graphs_action.setEnabled(False)

        self.open_temperature_graphs_action = QAction("Temperature Graphs", self)
        self.open_temperature_graphs_action.setStatusTip("Temperature Graphs")
        self.open_temperature_graphs_action.triggered.connect(self.open_temperature_graphs)
        self.open_temperature_graphs_action.setEnabled(False)

        self.open_scale_length_graphs_action = QAction("Gradient Scale Lengths", self)
        self.open_scale_length_graphs_action.setStatusTip("Gradient Scale Lengths")
        self.open_scale_length_graphs_action.triggered.connect(self.open_scale_length_graphs)
        self.open_scale_length_graphs_action.setEnabled(False)

        self.open_flux_graphs_action = QAction("Particle and Heat Fluxes", self)
        self.open_flux_graphs_action.setStatusTip("Particle and Heat Fluxes")
        self.open_flux_graphs_action.triggered.connect(self.open_flux_graphs)
        self.open_flux_graphs_action.setEnabled(False)

        self.open_electric_field_action = QAction("Electric field")
        self.open_electric_field_action.setStatusTip("Radial electric field")
        self.open_electric_field_action.triggered.connect(self.open_electric_field_graphs)
        self.open_electric_field_action.setEnabled(False)

        self.open_collisionality_action = QAction("Collisionality")
        self.open_collisionality_action.setStatusTip("Collisionality")
        self.open_collisionality_action.triggered.connect(self.open_collisionality_graphs)
        self.open_collisionality_action.setEnabled(False)

        self.open_density_graph_action = QAction("Densities")
        self.open_density_graph_action.setStatusTip("Density profiles")
        self.open_density_graph_action.triggered.connect(self.open_density_graphs)
        self.open_density_graph_action.setEnabled(False)


        exitAct = QAction(QIcon('exit.png'), '&Exit', self)
        exitAct.setShortcut('Ctrl+Q')
        exitAct.setStatusTip('Exit application')
        exitAct.triggered.connect(qApp.quit)

        #Add toolbar actions
        self.toolbar.addAction(self.open_GT3_action)
        self.toolbar.addAction(self.save_fig_action)

        menu = self.menuBar()
        menu.setNativeMenuBar(False)

        file_menu = menu.addMenu(u"&File")
        file_menu.addAction(self.open_GT3_action)
        file_menu.addAction(exitAct)

        edit_menu = menu.addMenu(u"&Edit")
        edit_menu.addAction(self.set_preferences_action)

        shot_menu = menu.addMenu(u"&Shot Data")
        shot_menu.addAction(self.open_velocity_graphs_action)
        shot_menu.addAction(self.open_temperature_graphs_action)
        shot_menu.addAction(self.open_scale_length_graphs_action)
        shot_menu.addAction(self.open_flux_graphs_action)
        shot_menu.addAction(self.open_electric_field_action)
        shot_menu.addAction(self.open_density_graph_action)
        shot_menu.addAction(self.open_collisionality_action)

    def open_velocity_graphs(self):

        velocities = [GraphContainer(self.chi_i.rho,
                                     [self.chi_i.shot.rtrans.vpol_D.val, self.chi_i.shot.rtrans.vpol_C.val],
                                     title="Poloidal Velocity",
                                     legend=[r"$V^D_{\theta}$", r"$V^C_{\theta}$"]),
                      GraphContainer(self.chi_i.rho,
                                     [self.chi_i.shot.rtrans.vtor_C_total.val, self.chi_i.shot.rtrans.vtor_D_total.val],
                                     title="Toroidal Velocity",
                                     legend=[r"$V^D_{\phi}$", r"$V^C_{\phi}$"])]
        velocity_graph = GraphPopout(velocities, grid=len(velocities), title="Shot Velocities")
        velocity_graph.show()

    def open_density_graphs(self):

        densities = [GraphContainer(self.chi_i.rho,
                                     [self.chi_i.shot.rtrans._n.D.val / (1E20), self.chi_i.shot.rtrans._n.e.val / (1E20)],
                                     title="Densities",
                                     legend=[r"$n_D$", r"$n_e$"])]
        density_graph = GraphPopout(densities, grid=len(densities), title="Density (x10^20 m^-3)")
        density_graph.show()

    def open_temperature_graphs(self):

        temperatures = [GraphContainer(self.chi_i.rho,
                                     [self.chi_i.shot.core.T.i.ev.fsa.val, self.chi_i.shot.core.T.e.ev.fsa.val],
                                     title="Temperatures",
                                     legend=[r"$T_i$", r"$T_e$"])]
        temperature_graph = GraphPopout(temperatures, grid=len(temperatures), title="Shot Temperatures")
        temperature_graph.show()

    def open_scale_length_graphs(self):

        scale_lengths = [GraphContainer(self.chi_i.rho,
                                        [self.chi_i.T.i.J.L.fsa.val, self.chi_i.T.e.J.L.fsa.val],
                                        title="Scale Lengths",
                                        legend=[r"$L_{T_i}$", r"$L_{T_e}$"]),
                         GraphContainer(self.chi_i.rho,
                                        [self.chi_i.n.i.L.fsa.val, self.chi_i.n.e.L.fsa.val],
                                        title="Scale Lengths",
                                        legend=[r"$L_{n_i}$", r"$L_{n_e}$"])
                         ]
        scale_length_graph = GraphPopout(scale_lengths, grid=len(scale_lengths), title="Gradient Scale Lengths")
        scale_length_graph.show()

    def open_flux_graphs(self):

        fluxes = [GraphContainer(self.chi_i.rho,
                                 [self.chi_i.shot.rtrans.gamma.D.diff.val / (1E20)],
                                 title="Ion Particle Flux",
                                 legend=[r"$\Gamma_{j,r}[10E20]$"]),
                  GraphContainer(self.chi_i.rho,
                                 [self.chi_i.shot.rtrans.Q.D.diff.val],
                                 title="Heat Flux",
                                 legend=[r"$Q_{r,j}$"])
                         ]
        flux_graphs = GraphPopout(fluxes, grid=len(fluxes), title="Particle/Heat fluxes")
        flux_graphs.show()

    def open_electric_field_graphs(self):

        field = [GraphContainer(self.chi_i.rho,
                                [self.chi_i.shot.core.E_r.fsa],
                                title="Electric Field",
                                legend=[r"$E_r$"])]
        field_graphs = GraphPopout(field, grid=len(field), title="Electric Field")
        field_graphs.show()

    def open_collisionality_graphs(self):

        field = [GraphContainer(self.chi_i.rho,
                                [self.chi_i.shot.rtrans.nustar],
                                title="Collisionality",
                                legend=[r"$\nu^*$"])]
        field_graphs = GraphPopout(field, grid=len(field), title="Collisionality")
        field_graphs.show()

    def saveFigureDialog(self):

        #f = QFileDialog.getOpenFileName(self, str("Save Figure"), os.getcwd(), str("Image Files (*)"))
        f = QFileDialog.getSaveFileName(self, str("Save Figure"), os.getcwd(), str("PNG (.png)"))
        self.mainGraph._saveFig(f[0])

    def set_preferences(self, triggered):
        settings = self.settings
        preferences_dialog = PreferencesDialog()
        preferences_dialog.setupUi(self)
        if preferences_dialog.exec_():
            print("YAY")
        else:
            print("FUCK")
        #default_config_value = settings.value(CONFIG_KEY_1, defaultValue=None, type=str)
        #preference_dialog = PreferencesDialog(default_config_value=default_config_value, parent=self)

if __name__=="__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
