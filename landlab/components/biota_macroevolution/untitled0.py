#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 13:36:34 2018

@author: njlyons
"""

fig, axes = plt.subplots(1, 1)

l=plt.Line2D([0,0.5],[1,1],color='k')
axes.add_line(l)

l.set_ydata([2,2])
#
#l.get_xdata()



#l=plt.Line2D([0,0.5],[2,2],linewidth=1,markeredgecolor='k',markerfacecolor='k',linestyle='-')
#axes.add_line(l)

axes.set_xlim([0,5])
axes.set_ylim([0,3])


n = '(((3:1000,4:1000)2:1000,(5:1000,6:1000)1:2000)0);'

for i, c in enumerate(n[::-1]):
    if c == ')':



