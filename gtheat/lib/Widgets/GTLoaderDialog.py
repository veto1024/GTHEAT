#!/usr/bin/pythnon

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from waitingspinnerwidget import QtWaitingSpinner

from ..chi_i import chi_i


class GT3loaderDialog(QWidget):
    def __init__(self, fileName, parent=None):
        super(GT3loaderDialog, self).__init__(parent)
        # self.initUI()
        self.spinner = QtWaitingSpinner(parent)

        self.spinner.setRoundness(70.0)
        self.spinner.setMinimumTrailOpacity(15.0)
        self.spinner.setTrailFadePercentage(70.0)
        self.spinner.setNumberOfLines(12)
        self.spinner.setLineLength(10)
        self.spinner.setLineWidth(5)
        self.spinner.setInnerRadius(10)
        self.spinner.setRevolutionsPerSecond(1)
        self.spinner.setColor(QColor(81, 4, 71))
        self.spinner.start()

        self.thread = QThread(self)
        self.worker = GT3Worker(fileName)
        self.worker.moveToThread(self.thread)  # worker will be run in another thread
        self.thread.started.connect(self.worker.doWork)  # Call worker.doWork when the thread starts

        self.worker.done.connect(self.stopSpinner)  # Stop the spinner once this is done.
        self.thread.start()  # Start the thread (and run doWork)
        self.thread.finished.emit()

    def stopSpinner(self):
        self.spinner.stop()
        self.shot = self.worker.chi_i


class GT3Worker(QObject):
    done = pyqtSignal()

    def __init__(self, fileName, parent=None):
        super(GT3Worker, self).__init__(parent)
        self.fileName = fileName

    def doWork(self):
        self.chi_i = chi_i(self.fileName[0])
        self.done.emit()