#!/usr/bin/python

from PyQt5 import QtCore, QtWidgets

class PreferencesDialog(QtWidgets.QDialog):

    def __init__(self, *args, **kwargs):
        super(PreferencesDialog, self).__init__(*args, **kwargs)

    def setupUi(self, parent):
        self.setParent(parent)
        self.setObjectName("Dialog")
        self.resize(700, 300)
        self._set_ok_reject_buttons()

        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)

        self.show()

    def _set_ok_reject_buttons(self):
        self.buttonBox = QtWidgets.QDialogButtonBox()
        self.buttonBox.setGeometry(QtCore.QRect(150, 250, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

