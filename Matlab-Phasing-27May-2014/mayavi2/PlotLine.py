"""

This module creates a 'Locator' axis, that can be used to mark a three
dimensional point in your data.

This code is distributed under the conditions of the BSD license.  See
LICENSE.txt for details.

Copyright (c) 2001-2002, Prabhu Ramachandran.
"""

__author__ = "Prabhu Ramachandran <prabhu_r@users.sf.net>"
__version__ = "$Revision: 1.5 $"
__date__ = "$Date: 2003/10/11 18:59:55 $"

import Base.Objects, Common
import Tkinter, tkColorChooser
import vtk as vtkpython
import vtkPipeline.vtkMethodParser
import math

debug = Common.debug

class PlotLine (Base.Objects.Module):

    """ This module creates a 'Locator' axis, that can be used to mark
    a three dimensional point in your data. """

    def __init__ (self, mod_m):
        debug ("In Locator::__init__ ()")
        Common.state.busy ()
        Base.Objects.Module.__init__ (self, mod_m)

        self.root = None
	self.axesactor = None
        self.slider = []
	self.resoln_var = []
	self.maxlen=0.0
        self.plotline_var = Tkinter.IntVar ()
        self.tubesize_var = Tkinter.DoubleVar ()
        self.linelenres_var = Tkinter.DoubleVar ()
	self.npoints_var = Tkinter.IntVar ()
        self.data_out = self.mod_m.GetOutput ()        

        self.line=vtkpython.vtkLineSource()
	self.probe=vtkpython.vtkProbeFilter()
        self.plotline = self.actor = vtkpython.vtkXYPlotActor ()

	self.tuber = vtkpython.vtkTubeFilter()
	self.lineactor = vtkpython.vtkActor()
	self.mapper = vtkpython.vtkPolyDataMapper()
        self.axes=vtkpython.vtkAxes()
	self.trans = vtkpython.vtkTransform()
        self.axesactor = vtkpython.vtkActor()

        self._initialize ()

        self.renwin.add_actors (self.lineactor)
        self.renwin.add_actors (self.plotline)
        self.pipe_objs = self.plotline
        self.renwin.Render ()
        Common.state.idle ()
        
    def _initialize (self):
        debug ("In Locator::_initialize ()")


        self.datacenter = []
        self.dim = []
	self.npoints_var.set(10)
	self.tubesize_var.set(10)
	self.linelenres_var.set(.1)
        self.coord = [0, 0, 0]
	self.curorien = [0.0,0.0,0.0]
        dim1 = self.data_out.GetBounds ()
	
        for i in range (3):
            self.dim.append ((dim1[2*i], dim1[2*i+1]))
            if dim1[2*i] != dim1[2*i+1]:
                self.coord[i] = 1
            val = dim1[2*i] + dim1[2*i+1]
            self.datacenter.append (val*0.5)

	self.maxlen = math.sqrt((self.dim[0][0]-self.dim[0][1])**2 +
				(self.dim[1][0]-self.dim[1][1])**2 +
			    	(self.dim[2][0]-self.dim[2][1])**2 )

	self.line.SetPoint1( (self.datacenter[0]-.25*self.maxlen, 
				self.datacenter[1],self.datacenter[2]) )
	self.line.SetPoint2( (self.datacenter[0]+.25*self.maxlen,
				self.datacenter[1],self.datacenter[2]) )

	self.line.SetResolution(self.npoints_var.get())
	self.set_tubesize()

	self.transline = vtkpython.vtkTransformPolyDataFilter()
	self.transline.SetInput( self.line.GetOutput() )
	self.transline.SetTransform( self.trans )

	self.probe.SetSource(self.data_out)
	self.probe.SetInput(self.transline.GetOutput())
	self.tuber.SetInput(self.probe.GetOutput())
	self.plotline.AddInput(self.probe.GetOutput())
	self.plotline.GetPositionCoordinate().SetValue(0.0, 0.0, 0)
	self.plotline.GetPosition2Coordinate().SetValue(1.0, 1.0, 0)
	self.plotline.SetXValuesToArcLength()
	self.plotline.GetProperty().SetColor(0,0,0)
	self.plotline.GetProperty().SetLineWidth(3)

        axmapper = vtkpython.vtkPolyDataMapper()
        axmapper.SetInput( self.axes.GetOutput() )

        self.axes.SymmetricOn()
        self.axesactor.SetMapper( axmapper )
        self.axesactor.GetProperty().SetLineWidth(self.tubesize_var.get()/2)
        self.axesactor.GetProperty ().SetColor (*Common.config.fg_color)
        self.axesactor.GetProperty().SetAmbient(1.0)
        self.axesactor.GetProperty().SetDiffuse(0.0)

        self.tuber.SetRadius( self.tubesize_var.get() )
        self.tuber.SetNumberOfSides( 12 )
	self.mapper.SetInput(self.tuber.GetOutput())
	self.lineactor.SetMapper(self.mapper)
        self.lineactor.GetProperty().SetAmbient(1.0)
        self.lineactor.GetProperty().SetDiffuse(0.0)


    def __del__ (self):
        debug ("In Locator::__del__ ()")
        if self.plotline:
            self.renwin.remove_actors (self.plotline)
            self.renwin.remove_actors (self.lineactor)
        self.renwin.Render ()        

    def SetInput (self, source):
        debug ("In Locator::SetInput ()")
        Common.state.busy ()
        self.data_out = source
        self._initialize ()
        Common.state.idle ()        

    def save_config (self, file): 
        debug ("In Locator::save_config ()")
        p = vtkPipeline.vtkMethodParser.VtkPickler ()
        for obj in (self.line, self.ax_map, self.plotline,
                    self.plotline.GetProperty ()):
            p.dump (obj, file)

    def load_config (self, file): 
        debug ("In Locator::load_config ()")
        p = vtkPipeline.vtkMethodParser.VtkPickler ()
        for obj in (self.axes, self.ax_map, self.plotline,
                    self.plotline.GetProperty ()):
            p.load (obj, file)

        self._initialize ()
        self.renwin.Render ()
        
    def config_changed (self): 
        debug ("In Locator::config_changed ()")
        self.plotline.GetProperty ().SetColor (*Common.config.fg_color)

    def make_main_gui (self):
        debug ("In Locator::make_main_gui ()")
        frame = Tkinter.Frame (self.root, relief='ridge', bd=2)
        frame.pack (side='top', pady=2, fill='both', expand=1)
        self.plotline_var.set (self.plotline.GetVisibility())
        plotline_but = Tkinter.Checkbutton (frame,
                                           text="Toggle Graph", 
                                           variable=self.plotline_var,
                                           onvalue=1, offvalue=0,
                                           command=self.set_plotline)
        plotline_but.grid (row=0, column=0, columnspan=3,
                          pady=4, sticky='w')

        rw = 1
        self.tubesize_var.set (self.tuber.GetRadius ())
        lab = Tkinter.Label (frame, text="Indicator Radius: ")
        lab.grid (row=rw, column=0, columnspan=1, sticky='w', pady=4)
        entr = Tkinter.Entry (frame, width=10, relief='sunken', 
                              textvariable=self.tubesize_var)
        entr.grid (row=rw, column=1, columnspan=2, sticky='w', pady=4)
        entr.bind ("<Return>", self.set_tubesize)
        rw = rw + 1

	lab = Tkinter.Label( frame, text="Sample Points: ")
	lab.grid(row=rw, column=0, columnspan=1, sticky='w', pady=4)
	entr = Tkinter.Entry( frame, width=10, relief='sunken',
				textvariable=self.npoints_var)
        entr.grid(row=rw, column=1, columnspan=2, sticky='w', pady=4)
	entr.bind ("<Return>", self.set_npoints)
	rw = rw + 1

#here
	p1=self.line.GetPoint1()
	p2=self.line.GetPoint2()
	curlinelen=abs(p1[0]-p2[0])
	sliderpos = curlinelen/self.maxlen
	self.linelen_slider = Tkinter.Scale(frame, label='Probe Length',
                              from_=0.0, to=1.0,
                              length="8c", orient='horizontal',
                              resolution=self.linelenres_var.get ())
	self.linelen_slider.set (sliderpos)
	self.linelen_slider.grid (row=rw, column=0, columnspan=3, sticky='ew', pady=4)
	self.linelen_slider.bind ("<ButtonRelease>", self.set_linelen)
	rw = rw + 1
	lab = Tkinter.Label( frame, text="Stretch Res: ")
	lab.grid(row=rw, column=0, columnspan=1, sticky='w', pady=4)
	entr = Tkinter.Entry( frame, width=10, relief='sunken',
				textvariable=self.linelenres_var)
        entr.grid(row=rw, column=1, columnspan=2, sticky='w', pady=4)
	entr.bind ("<Return>", self.set_linelenres)
	rw = rw + 1
 
	but = Tkinter.Button(frame, text="Move-Orient Line",
				underline=2, command=self.lineloc_gui)
	but.grid(row=rw, column=0, columnspan=2, sticky='ew')
	rw = rw + 1

        self.make_actor_gui (representation=0)


    def set_linelenres(self, event=None):
	self.linelen_slider.config(resolution=self.linelenres_var.get())

    def set_linelen(self, event=None):
	val=self.linelen_slider.get()
	self.line.SetPoint1( (self.datacenter[0] - .5*self.maxlen*val,
				self.datacenter[1], self.datacenter[2]) )
	self.line.SetPoint2( (self.datacenter[0] + .5*self.maxlen*val,
				self.datacenter[1], self.datacenter[2]) )
	self.renwin.Render()
	pass
	
    def lineloc_gui(self):

        sl_name = ["Point X", "Point Y", "Point Z", "Rotate X",
                    "Rotate Y", "Rotate Z"]

        self.ll_top = top = Tkinter.Toplevel (self.root)
        top.transient (self.root.master)
        top.protocol ("WM_DELETE_WINDOW", top.destroy)
        frame = Tkinter.Frame (top, relief='ridge', bd=2)
        frame.pack (side='top', fill='both', expand=1)
        rw = 0

	midpoint = self.lineactor.GetCenter()
	linexrange= self.lineactor.GetXRange()
	lineyrange= self.lineactor.GetYRange()
	linezrange= self.lineactor.GetZRange()


	for i in range(0,3):
	  if self.coord[i] == 1:
                self.resoln_var.append (Tkinter.DoubleVar ())
		sl = Tkinter.Scale(frame, label=sl_name[i],
                                    from_=self.dim[i][0], to=self.dim[i][1],
                                    length="8c", orient='horizontal',
                                    resolution=self.resoln_var[i].get ())
                sl.set (midpoint[i])
                sl.grid (row=rw, column=0, columnspan=3, sticky='ew')
                rw = rw + 1
                sl.bind ("<ButtonRelease>", self.translate_line)
                self.slider.append (sl)

	for i in range(3,6):
	  self.resoln_var.append(Tkinter.DoubleVar ())
          self.resoln_var[i].set (1.00)
	  if (self.coord[1] and self.coord[2]):
		sl =  Tkinter.Scale(frame, label=sl_name[i],
                                    from_=-180.0, to=180.0,
                                    length="8c", orient='horizontal',
                                    resolution=self.resoln_var[i].get ())
                sl.set (0.0)
                sl.grid (row=rw, column=0, columnspan=3, sticky='ew')
                rw = rw + 1
                sl.bind ("<ButtonRelease>", self.rotate_line)
                self.slider.append (sl)




	llen = 0.0
	for r in (linexrange, lineyrange,linezrange):
	  llen = llen + (r[0]-r[1])**2
	llen = math.sqrt( llen )
        self.axesactor.SetPosition( midpoint )
	self.axes.SetScaleFactor( llen*.25 )

	self.renwin.add_actors( self.axesactor )
	self.renwin.Render()

        but = Tkinter.Button (frame, text="Close", underline=0,
                              command=self.__lineloc_quit)
        but.grid(row=rw, column=1, sticky='ew')
        top.bind ("<Alt-c>", self.__lineloc_quit)

    def __lineloc_quit (self, event=None):
	self.renwin.remove_actors( self.axesactor )
	self.slider = []
	self.resoln_var = []
        self.ll_top.destroy ()
	self.renwin.Render()

    def translate_line(self, event=None):

	midpoint = self.lineactor.GetCenter()
        linexrange= self.lineactor.GetXRange()
        lineyrange= self.lineactor.GetYRange()
        linezrange= self.lineactor.GetZRange()

        delta = []
        for i in range (0, 3):
            if self.coord[i] == 1:
		sliderpos=self.slider[i].get()
		if (abs(sliderpos-midpoint[i]) > 0.0):
                  delta.append( sliderpos-midpoint[i] )
		else:
		  delta.append( 0.0 )
            else:
                delta.append(0.0)

	self.trans.Translate( *delta )

	midpoint = self.lineactor.GetCenter()
        self.axesactor.SetPosition(midpoint)
        self.renwin.Render ()

    def rotate_line(self, event=None):

	midpoint = self.lineactor.GetCenter()	
	minusmidpoint = [0,0,0]
	for i in range(3):
	  minusmidpoint[i] = -midpoint[i]

	delta = []
        for i in range (0, 3):
            if self.coord[i] == 1:
                sliderpos=self.slider[i+3].get()
                if (abs(sliderpos-self.curorien[i]) > 0.0):
                  delta.append( sliderpos-self.curorien[i] )
		  self.curorien[i] = sliderpos
                else:
                  delta.append( 0.0 )
            else:
                delta.append(0.0)

	self.trans.PostMultiply()
	self.trans.Translate( *minusmidpoint )
	self.transline.Update()

	self.trans.RotateX(delta[0])
	self.trans.RotateY(delta[1])
	self.trans.RotateZ(delta[2])
	self.trans.Translate( midpoint )
	self.transline.Update()

	midpoint = self.lineactor.GetCenter()	
        self.renwin.Render ()


    def set_plotline (self, event=None):
        debug ("In Locator::set_locator ()")
        self.plotline.SetVisibility (self.plotline_var.get ())
#        self.tuber.SetVisibility (self.plotline_var.get ())
        self.renwin.Render ()

    def set_npoints (self, event=None):
        debug ("In Locator::set_npoints ()")
	self.line.SetResolution(self.npoints_var.get())
        self.renwin.Render ()

    def set_tubesize (self, event=None):
        debug ("In Locator::set_tubesize ()")
	self.tuber.SetRadius(self.tubesize_var.get())
	if self.axesactor:
	  self.axesactor.GetProperty().SetLineWidth(self.tubesize_var.get()/5)
        self.renwin.Render ()

    def close_gui (self, event=None):
        debug ("In Locator::close_gui ()")
        Base.Objects.VizObject.close_gui (self, event)
        self.slider = []
