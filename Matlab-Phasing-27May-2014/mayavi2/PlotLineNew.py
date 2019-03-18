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

class PlotLineProbe:
    def __init__(self, renwin):

	self.renwin=renwin
	self.dim=[]
	self.locoriensliders=[]
	self.locorienslidervals=[0,0,0,0,0,0]
	self.locoriensliderres = [r0, r1, r2, r3, r4, r5] = \
		map(lambda resvars: Tkinter.DoubleVar(), range(6))
	self.coord=[0,0,0]

	self.line=vtkpython.vtkLineSource()
	self.tube=vtkpython.vtkTubeFilter()
	self.sphere1=vtkpython.vtkSphereSource()
	self.sphere2=vtkpython.vtkSphereSource()

	self.trans=vtkpython.vtkTransform()


	self.tline=vtkpython.vtkTransformPolyDataFilter()
	self.tline.SetTransform(self.trans)
	self.tline.SetInput(self.line.GetOutput())

	self.probe=vtkpython.vtkProbeFilter()
	self.probe.SetInput(self.tline.GetOutput())

	self.indicator=vtkpython.vtkAppendPolyData()
	self.tube.SetInput(self.probe.GetOutput())
	self.indicator.AddInput(self.tube.GetOutput())
	self.indicator.AddInput(self.sphere1.GetOutput())
	self.indicator.AddInput(self.sphere2.GetOutput())
	self.tindicator=vtkpython.vtkTransformPolyDataFilter()
	self.tindicator.SetInput(self.indicator.GetOutput())
	self.tindicator.SetTransform(self.trans)

	self.indicatoractor=vtkpython.vtkActor()
	self.indicatormapper=vtkpython.vtkPolyDataMapper()
	self.indicatormapper.SetInput(self.tindicator.GetOutput())
	self.indicatoractor.SetMapper(self.indicatormapper)

	pass

    def SetSource(self, source):
	self.data=source
	self.probe.SetSource(self.data)
	self.datacenter=[]
	dim1=self.data.GetBounds ()

        for i in range (3):
            self.dim.append ((dim1[2*i], dim1[2*i+1]))
            if dim1[2*i] != dim1[2*i+1]:
                self.coord[i] = 1
            val = dim1[2*i] + dim1[2*i+1]
            self.datacenter.append (val*0.5)

	self.locorienslidervals[3:6]=self.datacenter

        self.maxlen = math.sqrt((self.dim[0][0]-self.dim[0][1])**2 +
                                (self.dim[1][0]-self.dim[1][1])**2 +
                                (self.dim[2][0]-self.dim[2][1])**2 )

        self.line.SetPoint1( (self.datacenter[0]-.25*self.maxlen,
                                self.datacenter[1],self.datacenter[2]) )
	self.sphere1.SetCenter( (self.datacenter[0]-.25*self.maxlen,
                                self.datacenter[1],self.datacenter[2]) )
        self.line.SetPoint2( (self.datacenter[0]+.25*self.maxlen,
                                self.datacenter[1],self.datacenter[2]) )
	self.sphere2.SetCenter( (self.datacenter[0]+.25*self.maxlen,
                                self.datacenter[1],self.datacenter[2]) )

	pass

    def GetProbeOutput(self):
	return probe.GetOutput()
	pass

    def GetActor(self):
	return self.indicatoractor

    def config_gui(self, root=None):
        self.ll_top = top = Tkinter.Toplevel (root)
        frame = Tkinter.Frame (top, relief='ridge', bd=2)
        frame.pack (side='top', fill='both', expand=1)

	sl_name = ["Point X", "Point Y", "Point Z", "Rotate X",
                    "Rotate Y", "Rotate Z"]

        rw = 0
	for i in range(0,3):
		if self.coord[i]==1:
		  sl=Tkinter.Scale(frame, label=sl_name[i],
			from_=self.dim[i][0], to=self.dim[i][1],
			length="4c", orient='horizontal',
			resolution=self.locoriensliderres[i].get())
		  sl.set(self.locorienslidervals[i])
		  sl.grid(row=rw, column=0, columnspan=2, sticky='ew')	
		  sl.bind("<ButtonRelease>", self.translate)
		  self.locoriensliders.append(sl)
		  entr=Tkinter.Entry(frame, width=10,relief='sunken',
					textvariable=self.locoriensliderres[i])
		  entr.grid(row=rw, column=2, columnspan=1, sticky='ew')
		  entr.bind("<Return>", self.set_locoriensliderres(i))
		rw=rw+1

        but = Tkinter.Button (frame, text="Close", underline=0,
                               command=self.close_gui)
        but.grid(row=rw, column=0, columnspan=3, sticky='ew')
        top.bind ("<Alt-c>", self.close_gui)

    def close_gui(self):
	self.locoriensliders=[]
	self.ll_top.destroy()

    def translate(self, event=None):
        midpoint = self.indicatoractor.GetCenter()
        linexrange= self.indicatoractor.GetXRange()
        lineyrange= self.indicatoractor.GetYRange()
        linezrange= self.indicatoractor.GetZRange()

        delta = []
        for i in range (0, 3):
            if self.coord[i] == 1:
                sliderpos=self.locoriensliders[i].get()
                if (abs(sliderpos-midpoint[i]) > 0.0):
                  delta.append( sliderpos-midpoint[i] )
                else:
                  delta.append( 0.0 )
            else:
                delta.append(0.0)

        self.trans.Translate( *delta )

#        midpoint = self.indicatoractor.GetCenter()
#        self.axesactor.SetPosition(midpoint)
        self.renwin.Render ()

	pass

    def set_locoriensliderres(self, indx, event=None):
	self.locoriensliders[indx].config(resolution=\
					self.locoriensliderres[indx].get())
	

    def set_linelen(self, event=None):
        val=self.linelen_slider.get()
        self.line.SetPoint1( (self.datacenter[0] - .5*self.maxlen*val,
                                self.datacenter[1], self.datacenter[2]) )
        self.line.SetPoint2( (self.datacenter[0] + .5*self.maxlen*val,
                                self.datacenter[1], self.datacenter[2]) )
	self.renwin.Render()

class PlotLineNew (Base.Objects.Module):

    def __init__ (self, mod_m):
        debug ("In Locator::__init__ ()")
        Common.state.busy ()
        Base.Objects.Module.__init__ (self, mod_m)

        self.root = None
	self.axesactor = None
        self.slider = []
	self.resoln_var = []
	self.probes = []
	self.maxlen=0.0
        self.graph_var = Tkinter.IntVar ()
        self.tubesize_var = Tkinter.DoubleVar ()
        self.linelenres_var = Tkinter.DoubleVar ()
	self.npoints_var = Tkinter.IntVar ()
        self.data_out = self.mod_m.GetOutput ()        

        self.plot = self.actor = vtkpython.vtkXYPlotActor ()
	self.probes.append(PlotLineProbe(self.renwin))

        self._initialize ()

        self.renwin.add_actors (self.probes[0].GetActor())
#        self.renwin.add_actors (self.plot)
#        self.pipe_objs = self.plot
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

	self.probes[0].SetSource(self.data_out)

#	self.plot.AddInput(self.probes[0].GetLineOutput())
#	self.plot.GetPositionCoordinate().SetValue(0.0, 0.0, 0)
#	self.plot.GetPosition2Coordinate().SetValue(1.0, 1.0, 0)
#	self.plot.SetXValuesToArcLength()
#	self.plot.GetProperty().SetColor(0,0,0)
#	self.plot.GetProperty().SetLineWidth(3)

    def __del__ (self):
        debug ("In Locator::__del__ ()")
        if self.plot:
            self.renwin.remove_actors (self.plot)
	for p in self.probes:
            self.renwin.remove_actors (p.GetActor())
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
        for obj in (self.line, self.ax_map, self.plot,
                    self.plot.GetProperty ()):
            p.dump (obj, file)

    def load_config (self, file): 
        debug ("In Locator::load_config ()")
        p = vtkPipeline.vtkMethodParser.VtkPickler ()
        for obj in (self.axes, self.ax_map, self.plot,
                    self.plot.GetProperty ()):
            p.load (obj, file)

        self._initialize ()
        self.renwin.Render ()
        
    def config_changed (self): 
        debug ("In Locator::config_changed ()")
        self.plot.GetProperty ().SetColor (*Common.config.fg_color)

    def make_main_gui (self):
        frame = Tkinter.Frame (self.root, relief='ridge', bd=2)
        frame.pack (side='top', pady=2, fill='both', expand=1)

	rw=0
	but=Tkinter.Button(frame, text="Add Probe To Plot", 
						command=self.addplotlineprobe)
	but.grid(row=rw, column=0, columnspan=2, sticky='ew')
	rw=rw+1
	
	for p in self.probes:
	  pcbut=Tkinter.Button(frame, text="Config Probe"+str(rw-1),
							command=p.config_gui)
	  gcbut=Tkinter.Button(frame, text="Config Graph"+str(rw-1),
							command=self.nothing)
	  pcbut.grid(row=rw, column=0, columnspan=1, sticky='ew')
	  gcbut.grid(row=rw, column=1, columnspan=1, sticky='ew')
	  rw = rw + 1

    def nothing(self):
	print "nothing"
	pass

    def addplotlineprobe(self):
	self.probes.append(PlotLineProbe(self.renwin))
	self.close_gui()
	self.configure()
	pass
	
    def graph_gui(self):
	pass

    def set_plotline (self, event=None):
        debug ("In Locator::set_locator ()")
#        self.plot.SetVisibility (self.graph_var.get ())
#        self.tuber.SetVisibility (self.graph_var.get ())
        self.renwin.Render ()

    def set_npoints (self, event=None):
        debug ("In Locator::set_npoints ()")
	self.line.SetResolution(self.npoints_var.get())
        self.renwin.Render ()

    def close_gui (self, event=None):
        debug ("In Locator::close_gui ()")
        Base.Objects.VizObject.close_gui (self, event)

