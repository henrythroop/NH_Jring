#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 00:16:50 2018

This is a test to see if I can figure out how to do a simple GUI in Qt.

A: I'm close, but I get a lot of crashes. It's a lot faster than Tk to start up. But I'll probably just use Tk.

@author: throop
"""

import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QWidget, QLCDNumber, QSlider)

import random

import hbt
import os

class Window(QDialog):
#    https://stackoverflow.com/questions/12459811/how-to-embed-matplotlib-in-pyqt-for-dummies
    
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.button1 = QPushButton('Plot')
        self.button1.clicked.connect(self.plot)

        # Just some button connected to `imshow` method
        self.button2 = QPushButton('Imshow')
        self.button2.clicked.connect(self.imshow)

        # set the layout
        layout = QVBoxLayout()
#        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button1)
        layout.addWidget(self.button2)
        self.setLayout(layout)
        
        self.num = 10

    def plot(self):
        ''' plot some random stuff '''
        
        print("plot...")
        
        # random data
        data = [random.random() for i in range(self.num)]

        # instead of ax.hold(False)
        self.figure.clear()

        # create an axis
        ax = self.figure.add_subplot(111)

        # discards the old graph
        # ax.hold(False) # deprecated, see above

        # plot data
#        ax.imshow(hbt.dist(self.num))
        ax.plot(data, '*-')

        # refresh canvas
        self.canvas.draw()

    def imshow(self):
        ''' plot some random stuff '''

        print('Imshow...')
        # instead of ax.hold(False)
        self.figure.clear()

        # create an axis
        ax = self.figure.add_subplot(111)

        # plot data
        ax.imshow(hbt.dist(self.num))
#        ax.plot(data, '*-')

        # refresh canvas
        self.canvas.draw()

    def keyPressEvent(self, event):
        
        # http://zetcode.com/gui/pyqt5/eventssignals/
      print ("Key hit: {}".format(event.key()))
      
      if event.key() == Qt.Key_Escape:
         self.close()
      if event.key() == Qt.Key_Left:
          self.num -= 1
          self.plot()
      if event.key() == Qt.Key_Right:
          self.num += 1
          self.plot()

if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    main.show()

    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')

    sys.exit(app.exec_())