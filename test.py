# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 23:09:39 2016

@author: throop
"""

from Tkinter import *
from ttk import *

choices = ['1st Choice', '2nd Choice', '3rd Choice', '4th Choice']

root = Tk()

#for each in range(10):
OptionMenu(root, StringVar(), choices[0], *choices).grid()

root.mainloop()