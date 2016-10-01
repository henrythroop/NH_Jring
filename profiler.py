# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 00:12:37 2016

@author: throop
"""

import pstats

# Put this in my code:
#     cProfile.run('root.mainloop()', filename='profile.bin')  # This will call the __init__ function

# Then use this code to read and interpret that file.
# See https://docs.python.org/2/library/profile.html for more info.

p = pstats.Stats('profile.bin')

p.strip_dirs().sort_stats(2).print_stats(50)

p.sort_stats('cumulative').print_stats(50)
