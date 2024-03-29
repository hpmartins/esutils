#!/usr/bin/env python

# embedding_in_qt5.py --- Simple Qt5 application embedding matplotlib canvases
#
# Copyright (C) 2005 Florent Rougon
#			   2006 Darren Dale
#			   2015 Jens H Nielsen
#
# This file is an example program for matplotlib. It may be used and
# modified with no restriction; raw copies as well as modified versions
# may be distributed without limitation.

from __future__ import unicode_literals
import sys
import os
import random
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from numpy import arange, sin, pi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

progname = 'Singlecluster'
progversion = "0.1"


class MyMplCanvas(FigureCanvas):
	"""Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

	def __init__(self, parent=None, width=5, height=4, dpi=100):
		fig = Figure(figsize=(width, height), dpi=dpi)
		self.axes = fig.add_subplot(111)
		# We want the axes cleared every time plot() is called
		self.axes.hold(False)

		self.compute_initial_figure()

		#
		FigureCanvas.__init__(self, fig)
		self.setParent(parent)

		FigureCanvas.setSizePolicy(self,
								   QtWidgets.QSizePolicy.Expanding,
								   QtWidgets.QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)

	def compute_initial_figure(self):
		pass


class MyStaticMplCanvas(MyMplCanvas):
	"""Simple canvas with a sine plot."""

	def compute_initial_figure(self):
		t = arange(0.0, 3.0, 0.01)
		s = sin(2*pi*t)
		self.axes.plot(t, s)


class MyDynamicMplCanvas(MyMplCanvas):
	"""A canvas that updates itself every second with a new plot."""

	def __init__(self, *args, **kwargs):
		MyMplCanvas.__init__(self, *args, **kwargs)
		timer = QtCore.QTimer(self)
		timer.timeout.connect(self.update_figure)
		timer.start(1000)

	def compute_initial_figure(self):
		self.axes.plot([0, 1, 2, 3], [1, 2, 0, 4], 'r')

	def update_figure(self):
		# Build a list of 4 random integers between 0 and 10 (both inclusive)
		l = [random.randint(0, 10) for i in range(4)]

		self.axes.plot([0, 1, 2, 3], l, 'r')
		self.draw()


class ApplicationWindow(QtWidgets.QMainWindow):
	def __init__(self):
		QtWidgets.QMainWindow.__init__(self)
		self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
		self.setWindowTitle("application main window")

		# Adiciona File
		self.file_menu = QtWidgets.QMenu('&File', self)
		self.file_menu.addAction('&Open', self.fileOpen, QtCore.Qt.CTRL + QtCore.Qt.Key_O)
		self.file_menu.addAction('&Save', self.fileSave, QtCore.Qt.CTRL + QtCore.Qt.Key_S)
		self.file_menu.addAction('&Quit', self.fileQuit, QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
		self.menuBar().addMenu(self.file_menu)
		
		# Adiciona Reload
		self.reload_menu = QtWidgets.QMenu('&Reload', self)
		self.reload_menu.addAction('&Reload', self.reloadReload, QtCore.Qt.CTRL + QtCore.Qt.Key_R)
		self.menuBar().addMenu(self.reload_menu)

		# Adiciona Help
		self.help_menu = QtWidgets.QMenu('&Help', self)
		self.help_menu.addAction('&About', self.about)
		self.menuBar().addMenu(self.help_menu)

		# Main widget
		self.main_widget = QtWidgets.QWidget(self)

		l = QtWidgets.QVBoxLayout(self.main_widget)
		sc = MyStaticMplCanvas(self.main_widget, width=5, height=4, dpi=100)
		dc = MyDynamicMplCanvas(self.main_widget, width=5, height=4, dpi=100)
		l.addWidget(sc)
		l.addWidget(dc)

		self.main_widget.setFocus()
		self.setCentralWidget(self.main_widget)

		self.statusBar().showMessage("Status bar", 2000)

	def fileOpen(self):
		pass
	
	def fileSave(self):
		pass

	def fileQuit(self):
		self.close()
		
	def reloadReload(self):
		pass

	def closeEvent(self, ce):
		self.fileQuit()

	def about(self):
		QtWidgets.QMessageBox.about(self, "About", """TODO""")


qApp = QtWidgets.QApplication(sys.argv)
aw = ApplicationWindow()
aw.setWindowTitle("%s" % progname)
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()