#!/usr/bin/python

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import GT3
from GT3.Psi import LOWER_XPT, UPPER_XPT
from gtheat.lib.chi_i import chi_i

class GT3loaderDialog(QWidget):
    def __init__(self, fileName, parent=None):
        super(GT3loaderDialog, self).__init__(parent)
        self.fileName = fileName
        # self.initUI()
        self.get_run_kwargs()
        self.on_preferences_loaded()

    def on_preferences_loaded(self):
        self.GT3thread = QThread(self)
        self.worker = GT3Worker(self.fileName, pref_data=self.pref_data)

        self.worker.doWork()
        self.shot = self.worker.shot
        self.rho = self.shot.rtrans.rhor
        self.chi_i = self.worker.chi_i
        self.worker.moveToThread(self.GT3thread)  # worker will be run in another GT3thread
        self.GT3thread.started.connect(self.worker.doWork)  # Call worker.doWork when the GT3thread starts

        self.GT3thread.finished.emit()
        self.GT3thread.quit()

    # def start(self):
    #     self.GT3thread.start()  # Start the GT3thread (and run doWork)
    #     self.GT3thread.finished.connect(self.quitThread)
    #     print("Waiting now...")
    #     self.GT3thread.wait()
    #     print("Done executing waiting")


    def get_run_kwargs(self):
        thread = QThread(self)
        worker = GT3KwargObject()
        worker.doWork()
        worker.moveToThread(thread)
        thread.started.connect(worker.doWork)
        thread.finished.emit()
        thread.quit()
        self.pref_data = worker.data



class GT3KwargDialog(QDialog):
    def __init__(self, parent=None):
        super(GT3KwargDialog, self).__init__(parent)
        self.setWindowTitle("Run Options")
        self.parent = parent
        self.setupUi()

    def setupUi(self, parent=None):
        self.setParent(parent)
        self.setObjectName("GT3KwargDialog")
        self.resize(700, 300)
        self._set_ok_reject_buttons()
        self._set_group_boxes()
        self._set_group_labels()
        self._set_xpt_options()
        self._set_splines_options()

        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(10, 10, 10, 10)
        self.layout.addWidget(self.checkboxGroupsWidget)
        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)
        self.show()

    def _set_group_boxes(self):
        self.checkboxGroupsWidget = QWidget()
        self.checkboxGroupsLayout = QGridLayout()
        self.checkboxGroupsLayout.setGeometry(QRect(150, 250, 341, 32))
        self.checkboxGroupsWidget.setLayout(self.checkboxGroupsLayout)

    def _set_group_labels(self):
        self.checkboxGroupsLayout.addWidget(QLabel("Psi Options", self), 1, 1)
        self.checkboxGroupsLayout.addWidget(QLabel("Radial Transport Options", self), 1, 2)
        self.checkboxGroupsLayout.addWidget(QLabel("Spline Overrides", self), 1, 3)

    def _set_xpt_options(self):
        self.xptGroup = QGroupBox()
        self.xptGroupLayout = QVBoxLayout()
        self.xptGroupLabel = QLabel("X-point Override")
        self.xptGroupDefault = QRadioButton("Default")
        self.xptGroupUpper = QRadioButton("Upper X-Pt")
        self.xptGroupLower = QRadioButton("Lower X-Pt")
        self.xptGroupLayout.addWidget(self.xptGroupLabel)
        self.xptGroupLayout.addWidget(self.xptGroupDefault)
        self.xptGroupLayout.addWidget(self.xptGroupUpper)
        self.xptGroupLayout.addWidget(self.xptGroupLower)

        self.xptGroupDefault.setChecked(True)
        self.xptGroup.setContentsMargins(5,5,5,5)
        self.xptGroup.setLayout(self.xptGroupLayout)

        self.checkboxGroupsLayout.addWidget(self.xptGroup, 2, 1)

    def _set_splines_options(self):
        self.splineGroup = QGroupBox()
        self.splineGroupLayout = QVBoxLayout()
        self.splineGroupTi = QCheckBox(r"T_i")
        self.splineGroupTe = QCheckBox(r"T_e")
        self.splineGroupni = QCheckBox(r"n_i")
        self.splineGroupne = QCheckBox(r"n_e")
        self.splineGroupLayout.addWidget(self.splineGroupTe)
        self.splineGroupLayout.addWidget(self.splineGroupTi)
        self.splineGroupLayout.addWidget(self.splineGroupne)
        self.splineGroupLayout.addWidget(self.splineGroupni)
        self.splineGroup.setContentsMargins(5,5,5,5)
        self.splineGroup.setLayout(self.splineGroupLayout)
        self.checkboxGroupsLayout.addWidget(self.splineGroup, 2, 3)


    def _set_ok_reject_buttons(self):
        self.buttonBox = QDialogButtonBox()
        self.buttonBox.setGeometry(QRect(150, 250, 341, 32))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

    def _set_prefs(self):

        self.data = {
            "psi_args": {},
            "rtrans_override": {
                'T_i': False,
                'T_e': False,
                'n_i': False,
                'n_e': False,
            }
        }

        if self.xptGroupLower.isChecked():
            self.data["psi_args"] = {
                'xpt_select': LOWER_XPT
            }
        elif self.xptGroupUpper.isChecked():
            self.data["psi_args"] = {
                'xpt_select': UPPER_XPT
            }

        if self.splineGroupTi.isChecked():
            self.data["rtrans_override"]['T_i'] = True

        if self.splineGroupTe.isChecked():
            self.data["rtrans_override"]['T_e'] = True

        if self.splineGroupni.isChecked():
            self.data["rtrans_override"]['n_i'] = True

        if self.splineGroupne.isChecked():
            self.data["rtrans_override"]['n_e'] = True

    def accept(self):
        self._set_prefs()
        self.done(0)


class GT3KwargObject(QObject):
    done = pyqtSignal()
    def __init__(self, parent=None, **kwargs):
        super(GT3KwargObject, self).__init__(parent)
        self.parent = parent
        self.dialog = GT3KwargDialog()


    def doWork(self, parent=None):
        self.dialog.exec_()
        self.data = self.dialog.data


class GT3Worker(QObject):
    done = pyqtSignal()

    def __init__(self, fileName, parent=None, pref_data=None):
        super(GT3Worker, self).__init__(parent)
        self.fileName = fileName
        self.pref_data = pref_data

    def doWork(self, **kwargs):
        try:
            self.shot = GT3.gt3(self.fileName, **self.pref_data)
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


