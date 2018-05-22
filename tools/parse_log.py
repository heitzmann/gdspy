#!/bin/env python

######################################################################
#                                                                    #
#  Copyright 2009-2018 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import sys
import numpy
from matplotlib import pyplot

pts = {}
segs = {}

f1, a1 = pyplot.subplots(1, 1)
notes1 = []

with open('log.txt', 'r') as fin:
    for head in fin:
        if head == 'All pts:\n':
            data = fin.readline().strip().split()
            while len(data) > 0:
                lbl = data[0][-7:-1]
                x, y = [float(i) for i in data[1].split(',')]
                pts[lbl] = numpy.array((x, y))
                p, = a1.plot(x, y, '+')
                #a = a1.annotate(lbl, xy=(x,y), xycoords='data', ha='center')
                #a.set_visible(False)
                #notes1.append((p,a))
                data = fin.readline().strip().split()
        elif head == 'All segs:\n':
            data = fin.readline().strip().split()
            while len(data) > 0:
                lbl = data[0][-7:-1]
                segs[lbl] = (data[4][-6:], data[6][-7:-1])
                p1 = pts[segs[lbl][0]]
                p2 = pts[segs[lbl][1]]
                p, = a1.plot([p1[0], p2[0]], [p1[1], p2[1]], ':')
                #a = a1.annotate(lbl, xy=0.5*(p1+p2), xycoords='data', ha='center')
                #a.set_visible(False)
                #notes1.append((p,a))
                data = fin.readline().strip().split()
        elif head == 'Resulting segs:\n':
            data = fin.readline().strip().split()
            while len(data) > 0:
                lbl = data[0][-7:-1]
                segs[lbl] = (data[4][-6:], data[6][-7:-1])
                p1 = pts[segs[lbl][0]]
                p2 = pts[segs[lbl][1]]
                p, = a1.plot([p1[0], p2[0]], [p1[1], p2[1]], '--')
                #a = a1.annotate(lbl, xy=0.5*(p1+p2), xycoords='data', ha='center')
                #a.set_visible(False)
                #notes1.append((p,a))
                data = fin.readline().strip().split()


def on_move1(event):
    redraw = False
    for obj, note in notes1:
        vis = (obj.contains(event)[0] == True)
        if vis != note.get_visible():
            redraw = True
            note.set_visible(vis)
    if redraw:
        pyplot.draw()


#f1.canvas.mpl_connect('motion_notify_event', on_move1)
a1.grid(False)

pyplot.show()
