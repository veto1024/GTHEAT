
#!/usr/bin/python

# My python thing
import sys
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import time
from lib.chi_i import chi_i
import pyqtgraph as pg
from lib.Widgets.Color import Color
from lib.Widgets.GTLoaderDialog import GT3loaderDialog
import GT3
import os

class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setWindowTitle("GTHEAT 0.0.1 - A Cool Program")
        self.resize(800, 600)



        self.setContentsMargins(3, 3, 3, 3)

        layout = QGridLayout()

        # Set the toolbar
        self.toolbar = QToolBar("The Main Toolbar")

        #Set icon size
        self.toolbar.setIconSize(QSize(16, 16))

        # Create the toolbar
        self.toolbar = QToolBar("GTHEAT Toolbar")

        # Set icon size
        self.toolbar.setIconSize(QSize(16, 16))
        self.addToolBar(self.toolbar)
        # Add a button
        icon = QIcon("gtheat/lib/icons/bug.png")
        open_GT3_action = QAction("Open GT3 File", self)
        open_GT3_action.setIcon(icon)
        open_GT3_action.setStatusTip("Open a GT3 File")
        open_GT3_action.triggered.connect(self.onOpenGT3Started)
        self.toolbar.addAction(open_GT3_action)

        menu = self.menuBar()
        menu.setNativeMenuBar(False)

        file_menu = menu.addMenu(u"&File")
        file_menu.addAction(open_GT3_action)

        self.graphWidget = pg.PlotWidget()
        self.graphWidget.setBackground('w')
        self.graphWidget.setXRange(0.85, 1.0)
        self.setCentralWidget(self.graphWidget)


    def onOpenGT3Started(self, s):
        fileName = QFileDialog.getOpenFileName(self,
                                               str("Open File"), os.getcwd(), str("GT3 files (*)"))

        GT3loaderDialog(fileName, parent=self)
        while True:
            try:
                self.shot = GT3loaderDialog.shot  # type: GT3.gt3
                break
            except:
                pass
        self.graphWidget.plot(self.shot.rtrans.rhor, self.shot.rtrans.chi.i.chi4)

app = QApplication(sys.argv)
window = MainWindow()
window.show()
app.exec_()