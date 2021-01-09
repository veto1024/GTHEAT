#!/usr/bin/pythnon

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import GT3

from .chi_i import chi_i

class GT3loaderDialog(QWidget):
    def __init__(self, fileName, parent=None):
        super(GT3loaderDialog, self).__init__(parent)
        # self.initUI()


        self.thread = QThread(self)
        self.worker = GT3Worker(fileName)
        self.worker.doWork()
        self.shot = self.worker.shot
        self.rho = self.shot.rtrans.rhor
        self.chi_i = self.worker.chi_i
        self.worker.moveToThread(self.thread)  # worker will be run in another thread
        self.thread.started.connect(self.worker.doWork)  # Call worker.doWork when the thread starts

        self.thread.finished.emit()
        self.thread.quit()


    # def start(self):
    #     self.thread.start()  # Start the thread (and run doWork)
    #     self.thread.finished.connect(self.quitThread)
    #     print("Waiting now...")
    #     self.thread.wait()
    #     print("Done executing waiting")



class GT3Worker(QObject):
    done = pyqtSignal()

    def __init__(self, fileName, parent=None):
        super(GT3Worker, self).__init__(parent)
        self.fileName = fileName

    def doWork(self):
        try:
            self.shot = GT3.gt3(self.fileName)
        except:
            print("Invalid input file")
            self.done.emit()
        try:
            self.shot.run_radial_transport()
            self.shot.rtrans
        except Exception as e:
            print("Radial Transport failed to load")
            print("Error code: " + str(e))
            self.done.emit()
        try:
            self.chi_i = chi_i(self.shot)
            self.done.emit()
        except Exception as e:
            print("Could not load ion chi class")
            print("Error code: " + str(e))
            self.done.emit()


