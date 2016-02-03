########################################################################
##                                                                    ##
##  Copyright 2009-2016 Lucas Heitzmann Gabrielli                     ##
##                                                                    ##
##  This file is part of gdspy.                                       ##
##                                                                    ##
##  gdspy is free software: you can redistribute it and/or modify it  ##
##  under the terms of the GNU General Public License as published    ##
##  by the Free Software Foundation, either version 3 of the          ##
##  License, or any later version.                                    ##
##                                                                    ##
##  gdspy is distributed in the hope that it will be useful, but      ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of        ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ##
##  GNU General Public License for more details.                      ##
##                                                                    ##
##  You should have received a copy of the GNU General Public         ##
##  License along with gdspy.  If not, see                            ##
##  <http://www.gnu.org/licenses/>.                                   ##
##                                                                    ##
########################################################################

from __future__ import absolute_import

import os
import colorsys
import numpy
import Tkinter, tkMessageBox, tkColorChooser, tkFileDialog

import gdspy

__doc__ = """
Classes and functions for the visualization of layouts cerated with the
gdspy Python module.
"""

_stipple = tuple('@' + os.path.join(os.path.dirname(gdspy.__file__), 'data', '{:02d}.xbm'.format(n)) for n in range(10))
_icon_up = '@' + os.path.join(os.path.dirname(gdspy.__file__), 'data', 'up.xbm')
_icon_down = '@' + os.path.join(os.path.dirname(gdspy.__file__), 'data', 'down.xbm')
_invisible = 9
_icon_outline = '@' + os.path.join(os.path.dirname(gdspy.__file__), 'data', 'outline.xbm')

class ColorDict(dict):
    def __init__(self, default):
        super(ColorDict, self).__init__()
        self.default = default

    def __missing__(self, key):
        if self.default is None:
            layer, datatype = key
            rgb = '#{0[0]:02x}{0[1]:02x}{0[2]:02x}'.format([int(255 * c + 0.5) for c in colorsys.hsv_to_rgb((layer % 3)
                / 3.0 + (layer % 6 // 3) / 6.0 + (layer // 6) / 11.0, 1 - ((layer + datatype) % 8) / 12.0, 1 - (datatype % 3) / 4.0)])
        else:
            rgb = self.default
        self[key] = rgb
        return rgb


class PatternDict(dict):
    def __init__(self, default):
        super(PatternDict, self).__init__()
        self.default = default

    def __missing__(self, key):
        if self.default is None:
            pat = (key[0] + key[1]) % 8
        else:
            pat = self.default
        self[key] = pat
        return pat


class LayoutViewer(Tkinter.Frame):
    """
    Provide a GUI where the layout can be viewed.

    The view can be scrolled vertically with the mouse wheel, and
    horizontally by holding the shift key and using the mouse wheel.
    Dragging the 2nd mouse button also scrolls the view, and if control is
    held down, it scrolls 10 times faster.

    You can zoom in or out using control plus the mouse wheel, or drag a
    rectangle on the window with the 1st mouse button to zoom into that
    area.

    A ruler is available by clicking the 1st mouse button anywhere on the
    view and moving the mouse around.  The distance is shown in the status
    area.

    Double-clicking on any polygon gives some information about it.

    Color and pattern for each layer/datatype specification can be changed
    by left and right clicking on the icon in the layer/datatype list. Left
    and right clicking the text label changes the visibility.

    Parameters
    ----------
    cells : Cell, string or array-like
        The array of cells to be included in the view. If ``None``, all
        cells listed in ``Cell.cell_dict`` are used.
    hidden_types : array-like
        The array of tuples (layer, datatype) to start in hidden state.
    depth : integer
        Initial depth of referenced cells to be displayed.
    color : dictionary
        Dictionary of colors for each tuple (layer, datatype).  The colors
        must be strings in the format ``#rrggbb``.  A value with key
        ``default`` will be used as default color.
    pattern : dictionary
        Dictionary of patterns for each tuple (layer, datatype).  The
        patterns must be integers between 0 and 9, inclusive.  A value with
        key ``default`` will be used as default pattern.
    background : string
        Canvas background color in the format ``#rrggbb``.
    width : integer
        Horizontal size of the viewer canvas.
    height : integer
        Vertical size of the viewer canvas.

    Examples
    --------
    White background, filled shapes:

    >>> gdspy.LayoutViewer(pattern={'default':8}, background='#FFFFFF')

    No filling, black color for layer 0, datatype 1, automatic for others:

    >>> gdspy.LayoutViewer(pattern={'default':9}, color={(0,1):'#000000'})
    """
    def __init__(self, cells=None, hidden_types=[], depth=0, color={}, pattern={}, background='#202020', width=800, height=600):
        Tkinter.Frame.__init__(self, None)

        self.current_cell = Tkinter.StringVar()
        if cells is None:
            self.cells = gdspy.Cell.cell_dict
            cell_names = list(gdspy.Cell.cell_dict.keys())
            self.cell_bb = dict([(s, None) for s in self.cells])
        else:
            if isinstance(cells, str) or isinstance(cells, gdspy.Cell):
                cells = (cells,)
            self.cells = {}
            cell_names = []
            self.cell_bb = {}
            for c in cells:
                cell = gdspy.Cell.cell_dict.get(c, c)
                self.cells[cell.name] = cell
                cell_names.append(cell.name)
                self.cell_bb[cell.name] = None
        self.current_cell.set(cell_names[0])

        self.depth = Tkinter.IntVar()
        self.depth.set(depth)

        self.hidden_layers = hidden_types

        self.color = ColorDict(color.get('default', None))
        self.color.update(color)

        self.pattern = PatternDict(pattern.get('default', None))
        self.pattern.update(pattern)

        ## Setup resizable window
        self.grid(sticky='nsew')
        top = self.winfo_toplevel()
        top.rowconfigure(0, weight=1)
        top.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)

        ## Setup canvas
        self.canvas = Tkinter.Canvas(self, width=width, height=height, xscrollincrement=0, yscrollincrement=0)
        self.canvas.grid(row=0, column=0, sticky='nsew')
        self.canvas.configure(bg=background)
        bg = [int(c, 16) for c in (background[1:3], background[3:5], background[5:7])]
        self.default_outline = '#{0[0]:02x}{0[1]:02x}{0[2]:02x}'.format([(0 if c > 127 else 255) for c in bg])
        self.default_grey = '#{0[0]:02x}{0[1]:02x}{0[2]:02x}'.format([((c + 256) // 2 if c < 85 else (c // 2 if c > 171 else (255 if c > 127 else 0))) for c in bg])

        ## Setup scrollbars
        self.xscroll = Tkinter.Scrollbar(self, orient=Tkinter.HORIZONTAL, command=self.canvas.xview)
        self.xscroll.grid(row=1, column=0, sticky='ew')
        self.yscroll = Tkinter.Scrollbar(self, orient=Tkinter.VERTICAL, command=self.canvas.yview)
        self.yscroll.grid(row=0, column=1, sticky='ns')
        self.canvas['xscrollcommand'] = self.xscroll.set
        self.canvas['yscrollcommand'] = self.yscroll.set

        ## Setup toolbar
        self.frame = Tkinter.Frame(self)
        self.frame.columnconfigure(6, weight=1)
        self.frame.grid(row=2, column=0, columnspan=2, padx=2, pady=2, sticky='ew')

        ## Setup buttons
        self.home = Tkinter.Button(self.frame, text="Extents", command=self._update_canvas)
        self.zoom_in = Tkinter.Button(self.frame, text="In", command=self._zoom_in)
        self.zoom_out = Tkinter.Button(self.frame, text="Out", command=self._zoom_out)
        self.cell_menu = Tkinter.OptionMenu(self.frame, self.current_cell, *cell_names)
        self.depth_spin = Tkinter.Spinbox(self.frame, textvariable=self.depth, command=self._update_depth, from_=-1, to=128, increment=1, justify=Tkinter.RIGHT, width=3)
        self.home.grid(row=0, column=0, sticky='w')
        self.zoom_in.grid(row=0, column=1, sticky='w')
        self.zoom_out.grid(row=0, column=2, sticky='w')
        self.cell_menu.grid(row=0, column=3, sticky='w')
        Tkinter.Label(self.frame, text='Ref level:').grid(row=0, column=4, sticky='w')
        self.depth_spin.grid(row=0, column=5, sticky='w')
        self.bind_all("<KeyPress-Home>", self._update_canvas)
        self.bind_all("<KeyPress-A>", self._update_canvas)
        self.bind_all("<KeyPress-a>", self._update_canvas)

        ## Setup coordinates box
        self.coords = Tkinter.Label(self.frame, text="0, 0")
        self.coords.grid(row=0, column=6, sticky='e')

        ## Layers
        self.l_canvas = Tkinter.Canvas(self)
        self.l_canvas.grid(row=0, column=2, rowspan=3, sticky='nsew')
        self.l_scroll = Tkinter.Scrollbar(self, orient=Tkinter.VERTICAL, command=self.l_canvas.yview)
        self.l_scroll.grid(row=0, column=3, rowspan=3, sticky='ns')
        self.l_canvas['yscrollcommand'] = self.l_scroll.set

        ## Change current cell
        self.current_cell.trace_variable('w', self._update_canvas)

        ## Update depth
        self.depth_spin.bind('<FocusOut>', self._update_depth)

        ## Drag-scroll: button 2
        self.canvas.bind('<Button-2>', lambda evt: self.canvas.scan_mark(evt.x, evt.y))
        self.canvas.bind('<Motion>', self._mouse_move)

        ## Y scroll: scroll wheel
        self.canvas.bind('<MouseWheel>', lambda evt: self.canvas.yview(Tkinter.SCROLL, 1 if evt.delta < 0 else -1, Tkinter.UNITS))
        self.canvas.bind('<Button-4>', lambda evt: self.canvas.yview(Tkinter.SCROLL, -1, Tkinter.UNITS))
        self.canvas.bind('<Button-5>', lambda evt: self.canvas.yview(Tkinter.SCROLL, 1, Tkinter.UNITS))

        self.l_canvas.bind('<MouseWheel>', lambda evt: self.l_canvas.yview(Tkinter.SCROLL, 1 if evt.delta < 0 else -1, Tkinter.UNITS))
        self.l_canvas.bind('<Button-4>', lambda evt: self.l_canvas.yview(Tkinter.SCROLL, -1, Tkinter.UNITS))
        self.l_canvas.bind('<Button-5>', lambda evt: self.l_canvas.yview(Tkinter.SCROLL, 1, Tkinter.UNITS))

        ## X scroll: shift + scroll wheel
        self.bind_all('<Shift-MouseWheel>', lambda evt: self.canvas.xview(Tkinter.SCROLL, 1 if evt.delta < 0 else -1, Tkinter.UNITS))
        self.canvas.bind('<Shift-Button-4>', lambda evt: self.canvas.xview(Tkinter.SCROLL, -1, Tkinter.UNITS))
        self.canvas.bind('<Shift-Button-5>', lambda evt: self.canvas.xview(Tkinter.SCROLL, 1, Tkinter.UNITS))

        ## Object properties: double button 1
        ## Zoom rectangle: drag button 1
        ## Measure tool: button 1 (click + click, no drag)
        self.canvas.bind('<Button-1>', self._zoom_rect_mark)
        self.canvas.bind('<ButtonRelease-1>', self._mouse_btn_1)
        self.canvas.bind('<Double-Button-1>', self._properties)

        ## Zoom: control + scroll wheel
        self.bind_all('<Control-MouseWheel>', self._zoom)
        self.canvas.bind('<Control-Button-4>', self._zoom)
        self.canvas.bind('<Control-Button-5>', self._zoom)

        ## Update the viewer
        self.shown_cell = None
        self.shown_depth = depth
        self.canvas_margins = None
        self._update_canvas()
        self.master.title('gdspy - Layout Viewer')
        self.mainloop()

    def _update_depth(self, *args):
        try:
            d = self.depth.get()
        except:
            self.depth.set(self.shown_depth)
            return
        if d != self.shown_depth:
            self.shown_cell = self.current_cell.get()
            self.shown_depth = d
            pol_dict = self.cells[self.current_cell.get()].get_polygons(by_spec=True, depth=self.shown_depth)
            self._update_data(pol_dict)

    def _update_canvas(self, *args):
        pol_dict = self.cells[self.current_cell.get()].get_polygons(by_spec=True, depth=self.shown_depth)
        if self.shown_cell is None:
            width = float(self.canvas.cget('width'))
            height = float(self.canvas.cget('height'))
        else:
            width = float(self.canvas.winfo_width()) - self.canvas_margins[0]
            height = float(self.canvas.winfo_height()) - self.canvas_margins[1]
        self.shown_cell = self.current_cell.get()
        if self.cell_bb[self.current_cell.get()] is None:
            bb = [1e300, 1e300, -1e300, -1e300]
            for pols in pol_dict.itervalues():
                for pol in pols:
                    bb[0] = min(bb[0], pol[:,0].min())
                    bb[1] = min(bb[1], -pol[:,1].max())
                    bb[2] = max(bb[2], pol[:,0].max())
                    bb[3] = max(bb[3], -pol[:,1].min())
            self.cell_bb[self.current_cell.get()] = tuple(bb)
        else:
            bb = list(self.cell_bb[self.current_cell.get()])
        if bb[2] < bb[0]:
            tkMessageBox.showwarning('Warning', 'The selected cell is empty.')
            bb = [-1, -1, 1, 1]
        self.scale = ((bb[3] - bb[1]) / height, (bb[2] - bb[0]) / width)
        if self.scale[0] > self.scale[1]:
            self.scale = self.scale[0] * 1.05
            add = (width * self.scale - bb[2] + bb[0]) * 0.5
            bb[0] -= add
            bb[2] += add
            add = (bb[3] - bb[1]) * 0.025
            bb[1] -= add
            bb[3] += add
        else:
            self.scale = self.scale[1] * 1.05
            add = (height * self.scale - bb[3] + bb[1]) * 0.5
            bb[1] -= add
            bb[3] += add
            add = (bb[2] - bb[0]) * 0.025
            bb[0] -= add
            bb[2] += add
        self._update_data(pol_dict)
        self.canvas.configure(scrollregion=tuple([x / self.scale for x in bb]))
        self.canvas.zoom_rect = None
        if self.canvas_margins is None:
            self.update()
            self.canvas_margins = (int(self.canvas.winfo_width()) - width, int(self.canvas.winfo_height()) - height)

    def _update_data(self, pol_dict):
        self.canvas.delete(Tkinter.ALL)
        self.l_canvas.delete(Tkinter.ALL)
        self.canvas.ruler = None
        self.canvas.x_rl = 0
        self.canvas.y_rl = 0
        lbl_dict = dict([((label.layer, label.texttype), label) for label in self.cells[self.current_cell.get()].labels])
        layers = list(set(pol_dict.keys() + lbl_dict.keys()))
        layers.sort(reverse=True, key=lambda i: (-1, -1) if not isinstance(i, tuple) else i)
        self.l_canvas_info = []
        pos = 0
        wid = None
        hei = None
        bg = self.canvas['bg']
        self.l_canvas.configure(bg=bg)
        for i in layers:
            if i in self.hidden_layers:
                state = 'hidden'
                fg = self.default_grey
            else:
                state = 'normal'
                fg = self.default_outline
            if isinstance(i, tuple):
                lbl = (Tkinter.Label(self, bitmap=_icon_outline if self.pattern[i]==_invisible else _stipple[self.pattern[i]], bd=0, fg=self.color[i], bg=bg, anchor='c'),
                       Tkinter.Label(self, text='{0[0]}/{0[1]}'.format(i), bd=0, fg=fg, bg=bg, height=1, anchor='c', padx=8),
                       Tkinter.Label(self, bitmap=_icon_up, bd=0, fg=self.default_outline, bg=bg, anchor='c'),
                       Tkinter.Label(self, bitmap=_icon_down, bd=0, fg=self.default_outline, bg=bg, anchor='c'))
                lbl[0].bind('<Button-1>', self._change_color(lbl[0], i))
                lbl[0].bind('<Button-2>', self._change_pattern(lbl[0], i))
                lbl[0].bind('<Button-3>', self._change_pattern(lbl[0], i))
                lbl[1].bind('<Button-1>', self._change_visibility(lbl[1], i))
                lbl[1].bind('<Button-2>', self._change_other_visibility(i))
                lbl[1].bind('<Button-3>', self._change_other_visibility(i))
                lbl[2].bind('<Button-1>', self._raise(i))
                lbl[3].bind('<Button-1>', self._lower(i))
                for l in lbl:
                    l.bind('<MouseWheel>', lambda evt: self.l_canvas.yview(Tkinter.SCROLL, 1 if evt.delta < 0 else -1, Tkinter.UNITS))
                    l.bind('<Button-4>', lambda evt: self.l_canvas.yview(Tkinter.SCROLL, -1, Tkinter.UNITS))
                    l.bind('<Button-5>', lambda evt: self.l_canvas.yview(Tkinter.SCROLL, 1, Tkinter.UNITS))
                if wid is None:
                    lbl[1].configure(text='255/255')
                    hei = max(lbl[0].winfo_reqheight(), lbl[1].winfo_reqheight())
                    wid = lbl[1].winfo_reqwidth()
                    lbl[1].configure(text='{0[0]}/{0[1]}'.format(i))
                ids = (self.l_canvas.create_window(0, pos, window=lbl[0], anchor='sw'),
                       self.l_canvas.create_window(hei, pos, window=lbl[1], anchor='sw'),
                       self.l_canvas.create_window(hei + wid, pos, window=lbl[2], anchor='sw'),
                       self.l_canvas.create_window(2 * hei + wid, pos, window=lbl[3], anchor='sw'))
                self.l_canvas_info.append((i, ids, lbl))
                pos -= hei
            if i in pol_dict:
                if not isinstance(i, tuple):
                    for pol in pol_dict[i]:
                        self.canvas.create_polygon(*list((numpy.array((1, -1)) * pol / self.scale).flatten()), fill='',
                                outline=self.default_outline, activeoutline=self.default_outline, activewidth=2,
                                tag=('L' + str(i), 'V' + str(pol.shape[0])), state=state, dash=(8,8))
                        self.canvas.create_text(pol[:,0].mean() / self.scale, pol[:,1].mean() / -self.scale,
                                text=i, anchor=Tkinter.CENTER, fill=self.default_outline, tag=('L' + str(i), 'TEXT'))
                else:
                    for pol in pol_dict[i]:
                        self.canvas.create_polygon(*list((numpy.array((1, -1)) * pol / self.scale).flatten()),
                                fill=self.color[i], stipple=_stipple[self.pattern[i]],
                                offset='{},{}'.format(*numpy.random.randint(16, size=2)), outline=self.color[i],
                                activeoutline=self.default_outline, activewidth=2, tag=('L' + str(i), 'V' +
                                    str(pol.shape[0])), state=state)
            if i in lbl_dict:
                label = lbl_dict[i]
                self.canvas.create_text(label.position[0] / self.scale, label.position[1] / -self.scale,
                        text=label.text, anchor=[Tkinter.NW, Tkinter.N, Tkinter.NE, None, Tkinter.W, Tkinter.CENTER,
                            Tkinter.E, None, Tkinter.SW, Tkinter.S, Tkinter.SE][label.anchor], fill=self.color[i],
                        activefill=self.default_outline, tag=('L' + str(i), 'TEXT'))
        if (wid is None) or (hei is None) or (pos is None):
            pos = -12
            hei = 12
            wid = 12
        self.l_canvas.configure(width=3*hei+wid, scrollregion=(0,pos,3*hei+wid,0), yscrollincrement=hei)

    def _change_color(self, lbl, layer):
        def func(*args):
            rgb, color = tkColorChooser.askcolor(self.color[layer], title='Select color')
            if not color is None:
                self.color[layer] = color
                lbl.configure(fg=color)
                for i in self.canvas.find_withtag('L' + str(layer)):
                    self.canvas.itemconfigure(i, fill=color)
                    if layer[0] >= 0 and 'TEXT' not in self.canvas.gettags(i):
                        self.canvas.itemconfigure(i, outline=color)
        return func

    def _change_pattern(self, lbl, layer):
        top = self.winfo_toplevel()
        def func(*args):
            pattern = []
            dlg = Tkinter.Toplevel()
            dlg.title('Select pattern')
            dlg.resizable(False, False)
            for i in range(10):
                choice = Tkinter.Button(dlg, bitmap=_stipple[i], command=(lambda x: (lambda: pattern.append(x) or
                    dlg.destroy()))(i))
                choice.grid(row=0, column=i, padx=3, pady=3)
            choice = Tkinter.Button(dlg, text='Cancel', command=dlg.destroy)
            choice.grid(row=1, column=0, columnspan=10, padx=3, pady=3, sticky='e')
            dlg.focus_set()
            dlg.wait_visibility()
            dlg.grab_set()
            dlg.wait_window(dlg)
            if len(pattern) > 0:
                self.pattern[layer] = pattern[0]
                lbl.configure(bitmap=_icon_outline if pattern[0]==_invisible else _stipple[pattern[0]])
                for i in self.canvas.find_withtag('L' + str(layer)):
                    if 'TEXT' not in self.canvas.gettags(i):
                        self.canvas.itemconfigure(i, stipple=_stipple[pattern[0]])
        return func

    def _change_visibility(self, lbl, layer):
        def func(*args):
            if layer in self.hidden_layers:
                self.hidden_layers.remove(layer)
                lbl.configure(fg=self.default_outline)
                for j in self.canvas.find_withtag('L' + str(layer)):
                    self.canvas.itemconfigure(j, state='normal')
            else:
                self.hidden_layers.append(layer)
                lbl.configure(fg=self.default_grey)
                for j in self.canvas.find_withtag('L' + str(layer)):
                    self.canvas.itemconfigure(j, state='hidden')
        return func

    def _change_other_visibility(self, layer):
        def func(*args):
            for other, ids, lbl in self.l_canvas_info:
                if layer != other:
                    unhide = other in self.hidden_layers
                    break
            for other, ids, lbl in self.l_canvas_info:
                if layer != other:
                    if unhide and (other in self.hidden_layers):
                        self.hidden_layers.remove(other)
                        lbl[1].configure(fg=self.default_outline)
                        for j in self.canvas.find_withtag('L' + str(other)):
                            self.canvas.itemconfigure(j, state='normal')
                    elif not (unhide or (other in self.hidden_layers)):
                        self.hidden_layers.append(other)
                        lbl[1].configure(fg=self.default_grey)
                        for j in self.canvas.find_withtag('L' + str(other)):
                            self.canvas.itemconfigure(j, state='hidden')
        return func

    def _raise(self, layer):
        def func(*args):
            idx = 0
            while self.l_canvas_info[idx][0] != layer:
                idx += 1
            if idx < len(self.l_canvas_info) - 1:
                hei = int(self.l_canvas['yscrollincrement'])
                under, idu, _ = self.l_canvas_info[idx]
                above, ida, _ = self.l_canvas_info[idx + 1]
                self.canvas.tag_raise('L' + str(under), 'L' + str(above))
                for i in idu:
                    self.l_canvas.move(i, 0, -hei)
                for i in ida:
                    self.l_canvas.move(i, 0, hei)
                self.l_canvas.yview_scroll(-1, Tkinter.UNITS)
                self.l_canvas_info[idx], self.l_canvas_info[idx + 1] = self.l_canvas_info[idx + 1], self.l_canvas_info[idx]
        return func

    def _lower(self, layer):
        def func(*args):
            idx = 0
            while self.l_canvas_info[idx][0] != layer:
                idx += 1
            if idx > 0:
                hei = int(self.l_canvas['yscrollincrement'])
                under, idu, _ = self.l_canvas_info[idx - 1]
                above, ida, _ = self.l_canvas_info[idx]
                self.canvas.tag_raise('L' + str(under), 'L' + str(above))
                for i in idu:
                    self.l_canvas.move(i, 0, -hei)
                for i in ida:
                    self.l_canvas.move(i, 0, hei)
                self.l_canvas.yview_scroll(1, Tkinter.UNITS)
                self.l_canvas_info[idx], self.l_canvas_info[idx - 1] = self.l_canvas_info[idx - 1], self.l_canvas_info[idx]
        return func

    def _mouse_move(self, evt):
        x = self.canvas.canvasx(evt.x)
        y = self.canvas.canvasy(evt.y)
        if self.canvas.ruler is None:
            self.coords.configure(text='{0:g}, {1:g}'.format(x * self.scale, -y * self.scale))
        else:
            self.canvas.coords(self.canvas.ruler, self.canvas.x_rl, self.canvas.y_rl, x, y)
            dx = (x - self.canvas.x_rl) * self.scale
            dy = (self.canvas.y_rl - y) * self.scale
            self.coords.configure(text='Distance: {0:g} | dx = {1:g} | dy = {2:g}'.format((dx**2 + dy**2)**0.5, dx, dy))
        if int(evt.state) & 0x0200:
            if int(evt.state) & 0x0004:
                self.canvas.scan_dragto(evt.x, evt.y, 10)
            else:
                self.canvas.scan_dragto(evt.x, evt.y, 1)
        elif int(evt.state) & 0x0100:
            if self.canvas.zoom_rect is None:
                self.canvas.zoom_rect = self.canvas.create_rectangle(self.canvas.x_zr, self.canvas.y_zr, self.canvas.x_zr, self.canvas.y_zr, outline='#DDD')
            self.canvas.coords(self.canvas.zoom_rect, self.canvas.x_zr, self.canvas.y_zr, x, y)

    def _zoom(self, evt):
        if evt.num==4:
            evt.delta = 1
        elif evt.num==5:
            evt.delta = -1
        s = 1.5 if evt.delta > 0 else 1/1.5
        self.scale /= s
        x0 = s * self.canvas.canvasx(evt.x) - evt.x
        y0 = s * self.canvas.canvasy(evt.y) - evt.y
        self.canvas.scale(Tkinter.ALL, 0, 0, s, s)
        self.canvas.x_rl *= s
        self.canvas.y_rl *= s
        bb = self.canvas.bbox(Tkinter.ALL)
        if not bb is None:
            w = (bb[2] - bb[0]) * 1.2
            h = (bb[3] - bb[1]) * 1.2
            bb = (bb[0] - w, bb[1] - h, bb[2] + w, bb[3] + h)
            self.canvas['scrollregion'] = bb
            self.canvas.xview(Tkinter.MOVETO, (x0 - bb[0]) / (bb[2] - bb[0]))
            self.canvas.yview(Tkinter.MOVETO, (y0 - bb[1]) / (bb[3] - bb[1]))

    def _zoom_in(self):
        s = 1.5
        self.scale /= s
        self.canvas.scale(Tkinter.ALL, 0, 0, s, s)
        self.canvas.x_rl *= s
        self.canvas.y_rl *= s
        bb = self.canvas.bbox(Tkinter.ALL)
        w = (bb[2] - bb[0]) * 1.2
        h = (bb[3] - bb[1]) * 1.2
        bb = (bb[0] - w, bb[1] - h, bb[2] + w, bb[3] + h)
        self.canvas['scrollregion'] = bb
        x0 = self.xscroll.get()
        x0 = 0.5 * (x0[1] + x0[0]) - 0.5 * (float(self.canvas.winfo_width()) - self.canvas_margins[0]) / (bb[2] - bb[0])
        y0 = self.yscroll.get()
        y0 = 0.5 * (y0[1] + y0[0]) - 0.5 * (float(self.canvas.winfo_height()) - self.canvas_margins[1]) / (bb[3] - bb[1])
        self.canvas.xview(Tkinter.MOVETO, x0)
        self.canvas.yview(Tkinter.MOVETO, y0)

    def _zoom_out(self):
        s = 1/1.5
        self.scale /= s
        self.canvas.scale(Tkinter.ALL, 0, 0, s, s)
        self.canvas.x_rl *= s
        self.canvas.y_rl *= s
        bb = self.canvas.bbox(Tkinter.ALL)
        w = (bb[2] - bb[0]) * 1.2
        h = (bb[3] - bb[1]) * 1.2
        bb = (bb[0] - w, bb[1] - h, bb[2] + w, bb[3] + h)
        self.canvas['scrollregion'] = bb
        x0 = self.xscroll.get()
        x0 = 0.5 * (x0[1] + x0[0]) - 0.5 * (float(self.canvas.winfo_width()) - self.canvas_margins[0]) / (bb[2] - bb[0])
        y0 = self.yscroll.get()
        y0 = 0.5 * (y0[1] + y0[0]) - 0.5 * (float(self.canvas.winfo_height()) - self.canvas_margins[1]) / (bb[3] - bb[1])
        self.canvas.xview(Tkinter.MOVETO, x0)
        self.canvas.yview(Tkinter.MOVETO, y0)

    def _zoom_rect_mark(self, evt):
        self.canvas.x_zr = float(self.canvas.canvasx(evt.x))
        self.canvas.y_zr = float(self.canvas.canvasy(evt.y))

    def _mouse_btn_1(self, evt):
        if self.canvas.zoom_rect is None:
            if self.canvas.ruler is None:
                x0 = self.canvas.canvasx(evt.x)
                y0 = self.canvas.canvasy(evt.y)
                self.canvas.ruler = self.canvas.create_line(x0, y0, x0, y0, arrow=Tkinter.BOTH, fill=self.default_outline, width=2)
                self.canvas.x_rl = x0
                self.canvas.y_rl = y0
            else:
                self.canvas.delete(self.canvas.ruler)
                self.canvas.ruler = None
        else:
            x1 = float(self.canvas.winfo_width()) - self.canvas_margins[0]
            sx = float(self.canvas.canvasx(evt.x))
            dx = abs(self.canvas.x_zr - sx)
            sx += self.canvas.x_zr
            y1 = float(self.canvas.winfo_height()) - self.canvas_margins[1]
            sy = float(self.canvas.canvasy(evt.y))
            dy = abs(self.canvas.y_zr - sy)
            sy += self.canvas.y_zr
            self.canvas.delete(self.canvas.zoom_rect)
            self.canvas.zoom_rect = None
            if abs(dx * dy) > 1.0e-12:
                s = (x1 / dx, y1 / dy)
                if s[0] < s[1]:
                    s = s[0]
                    y0 = 0.5 * (s * sy - y1)
                    x0 = 0.5 * s * (sx - dx)
                else:
                    s = s[1]
                    x0 = 0.5 * (s * sx - x1)
                    y0 = 0.5 * s * (sy - dy)
                self.scale /= s
                self.canvas.scale(Tkinter.ALL, 0, 0, s, s)
                self.canvas.x_rl *= s
                self.canvas.y_rl *= s
                bb = self.canvas.bbox(Tkinter.ALL)
                if not bb is None:
                    w = (bb[2] - bb[0]) * 1.5
                    h = (bb[3] - bb[1]) * 1.5
                    bb = (bb[0] - w, bb[1] - h, bb[2] + w, bb[3] + h)
                    self.canvas['scrollregion'] = bb
                    self.canvas.xview(Tkinter.MOVETO, (x0 - bb[0]) / (bb[2] - bb[0]))
                    self.canvas.yview(Tkinter.MOVETO, (y0 - bb[1]) / (bb[3] - bb[1]))

    def _properties(self, evt):
        if self.canvas.ruler is not None:
            self.canvas.delete(self.canvas.ruler)
            self.canvas.ruler = -1
        i = self.canvas.find_closest(self.canvas.canvasx(evt.x), self.canvas.canvasy(evt.y))
        bb = self.canvas.bbox(i)
        if not bb is None:
            bb = (bb[0] * self.scale, -bb[3] * self.scale, bb[2] * self.scale, -bb[1] * self.scale)
            tags = self.canvas.gettags(i)
            if not 'TEXT' in tags:
                tkMessageBox.showinfo('Element information', 'Layer/datatpe: {0}\nVertices: {1}\nApproximate bounding box:\n({2[0]:g}, {2[1]:g}) - ({2[2]:g}, {2[3]:g})'.format(tags[0][1:], tags[1][1:], bb), parent=self.canvas)
