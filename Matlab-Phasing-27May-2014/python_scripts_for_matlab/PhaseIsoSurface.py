"""

This module shows an iso-surface of scalar data. This will work for
any dataset.

This code is distributed under the conditions of the BSD license.  See
LICENSE.txt for details.

Copyright (c) 2001-2002, Prabhu Ramachandran.
"""

__author__ = "Prabhu Ramachandran <prabhu_r@users.sf.net>"
__version__ = "$Revision: 1.4 $"
__date__ = "$Date: 2004/07/06 04:29:31 $"

import Base.Objects, Common
import Tkinter, tkColorChooser, math
import vtk as vtkpython
import vtkPipeline.vtkMethodParser

debug = Common.debug

class PhaseIsoSurface (Base.Objects.Module):

    """ This module shows an iso-surface of scalar data. This will
    work for any dataset."""

    def __init__ (self, mod_m): 
        debug ("In IsoSurface::__init__ ()")
        Common.state.busy ()
        Base.Objects.Module.__init__ (self, mod_m)
        self.cont_fil = vtkpython.vtkContourFilter ()
        self.norm = vtkpython.vtkPolyDataNormals ()
        self.mapper = self.map = vtkpython.vtkPolyDataMapper ()
        self.actor = self.act = vtkpython.vtkActor ()
        self.data_out = self.mod_m.GetOutput ()
        self._initialize ()
        self._gui_init ()
        if not self.mod_m.get_scalar_data_name ():
            msg = "Warning: No scalar data present to contour!"
            Common.print_err (msg)
        self.renwin.Render ()
        Common.state.idle ()

    def __del__ (self): 
        debug ("In IsoSurface::__del__ ()")
        if self.act:
            self.renwin.remove_actors (self.act)
        self.renwin.Render ()

    def _initialize (self):
        debug ("In IsoSurface::_initialize ()")
        self.cont_fil.SetInput (self.data_out)

        self.cont_fil.ComputeScalarsOff()
        self.cont_fil.Update()
        isopolydata=self.cont_fil.GetOutput()
        isopolydata.Update()
        print isopolydata.GetPointData()

        self.map.SetInput (self.cont_fil.GetOutput ())
        self.map.SetLookupTable (self.mod_m.get_scalar_lut ())
        dr = self.mod_m.get_scalar_data_range ()
        self.data_range = dr
#        self.map.SetScalarRange (dr)

	self.map.SetScalarModeToUsePointFieldData()
	self.map.SetScalarRange( (-3.1415,3.1415) )
	self.map.ColorByArrayComponent( "phases", 0)

        self.act.SetMapper (self.map)
        #self.act.GetProperty ().BackfaceCullingOff ()
        #self.act.GetProperty ().FrontfaceCullingOff ()
        midrange = (dr[0] + dr[1])*0.5
        self.slider_pos = midrange
        self.cont_fil.SetValue (0, midrange)
        self.renwin.add_actors (self.act)
        # used for the pipeline browser
        self.pipe_objs = self.act
        
    def _gui_init (self): 
        debug ("In IsoSurface::_gui_init ()")
        self.slider = None
        self.contour_var = Tkinter.DoubleVar ()
        self.contour_var.set (self.slider_pos)
        self.c_entry_var = Tkinter.DoubleVar ()
        self.c_entry_var.set (self.slider_pos)
        dr = self.mod_m.get_scalar_data_range ()
        self.min_cnt = Tkinter.DoubleVar ()
        self.max_cnt = Tkinter.DoubleVar ()
        self.min_cnt.set (dr[0])
        self.max_cnt.set (dr[1])
	self.min_phcolor = Tkinter.DoubleVar ()
	self.max_phcolor = Tkinter.DoubleVar ()
	self.min_phcolor.set(-3.1415)
	self.max_phcolor.set(3.1415)
        self.normals_on = Tkinter.IntVar ()
        self.normals_on.set (0)
        self.angle_var = Tkinter.DoubleVar ()
        self.angle_var.set (45)
        self._auto_sweep_init ()
        self.sweep_step.set (15)

    def SetInput (self, source): 
        debug ("In IsoSurface::SetInput ()")
        Common.state.busy ()
        self.data_out = source
        self.cont_fil.SetInput (self.data_out)
        dr = self.mod_m.get_scalar_data_range ()
        if (dr[0] != self.data_range[0]) or (dr[1] != self.data_range[1]):
            self.data_range = dr
            self.map.SetScalarRange (dr)
            self.min_cnt.set (dr[0])
            self.max_cnt.set (dr[1])        
            self.c_entry_var.set ((dr[0] + dr[1])*0.5)
            self.change_contour_limits ()
            self.change_contour_entry ()
        Common.state.idle ()

    def save_config (self, file): 
        debug ("In IsoSurface::save_config ()")
        file.write ("%f, %d, %f, %f, %f\n"%(self.c_entry_var.get (),
                                            self.normals_on.get (),
                                            self.angle_var.get (),
                                            self.min_cnt.get (),
                                            self.max_cnt.get ()))
        p = vtkPipeline.vtkMethodParser.VtkPickler ()
        for obj in (self.cont_fil, self.norm, self.map,
                    self.act, self.act.GetProperty ()):
            p.dump (obj, file)

    def load_config (self, file): 
        debug ("In IsoSurface::load_config ()")
        c_val, norm_on, ang, min_cnt, max_cnt = eval (file.readline ())
        self.slider_pos  = c_val
        self.contour_var.set (c_val)
        self.c_entry_var.set (c_val)
        self.normals_on.set (norm_on)
        self.angle_var.set (ang)
        self.min_cnt.set (min_cnt)
        self.max_cnt.set (max_cnt)
        p = vtkPipeline.vtkMethodParser.VtkPickler ()
        for obj in (self.cont_fil, self.norm, self.map,
                    self.act, self.act.GetProperty ()):        
            p.load (obj, file)

        self.change_contour_slider ()
        self.do_normals ()
        
    def config_changed (self): 
        debug ("In IsoSurface::config_changed ()")
        pass

    def make_main_gui (self): 
        debug ("In IsoSurface::make_main_gui ()")
        frame = Tkinter.Frame (self.root, relief='ridge', bd=2)
        frame.pack (side='top', fill='both', expand=1)

        name = "Scalar Variable: " + self.mod_m.get_scalar_data_name ()
        dr = (self.min_cnt.get (), self.max_cnt.get ())
        resolution = (dr[1] - dr[0])/1000
        self.slider = Tkinter.Scale (frame, label=name,
                                     resolution=resolution,
                                     variable=self.contour_var,
                                     from_=dr[0], to=dr[1], length="8c",
                                     orient='horizontal')

        self.slider.grid (row=0, column=0, columnspan=2, sticky='ew')
        self.slider.bind ("<ButtonRelease>", self.change_contour_slider )

        lab = Tkinter.Label (frame, text="Iso-surface value :")
        lab.grid (row=1, column=0, sticky='w')
        entry = Tkinter.Entry (frame, width=10, relief='sunken',
                               textvariable=self.c_entry_var)
        entry.grid (row=1, column=1, sticky='ew')
        entry.bind ("<Return>", self.change_contour_entry)

        lab = Tkinter.Label (frame, text="Minimum contour:")
        lab.grid (row=2, column=0, sticky='w')
        entry = Tkinter.Entry (frame, width=10, relief='sunken', 
                               textvariable=self.min_cnt)
        entry.grid (row=2, column=1, sticky='we')
        entry.bind ("<Return>", self.change_contour_limits)

        lab = Tkinter.Label (frame, text="Maximum contour:")
        lab.grid (row=3, column=0, sticky='w')
        entry = Tkinter.Entry (frame, width=10, relief='sunken', 
                               textvariable=self.max_cnt)
        entry.grid (row=3, column=1, sticky='we')
        entry.bind ("<Return>", self.change_contour_limits)

        lab = Tkinter.Label (frame, text="Minimum Phase:")
        lab.grid (row=4, column=0, sticky='w')
        entry = Tkinter.Entry (frame, width=10, relief='sunken', 
                               textvariable=self.min_phcolor)
        entry.grid (row=4, column=1, sticky='we')
        entry.bind ("<Return>", self.change_phase_scalarlimits)

        lab = Tkinter.Label (frame, text="Maximum Phase:")
        lab.grid (row=5, column=0, sticky='w')
        entry = Tkinter.Entry (frame, width=10, relief='sunken', 
                               textvariable=self.max_phcolor)
        entry.grid (row=5, column=1, sticky='we')
        entry.bind ("<Return>", self.change_phase_scalarlimits)

        norm = Tkinter.Checkbutton (frame, text="Compute PolyDataNormals", 
                                    variable=self.normals_on, onvalue=1,
                                    offvalue=0, command=self.do_normals)
        norm.grid (row=6, columnspan=2, sticky='w')
        lab = Tkinter.Label (frame, text="Feature Angle (in degrees): ")
        lab.grid (row=7, column=0, sticky='w')
        entry = Tkinter.Entry (frame, width=10, relief='sunken', 
                               textvariable=self.angle_var)
        entry.grid (row=7, column=1, sticky='ew')
        entry.bind ("<Return>", self.do_normals)

        self.make_actor_gui ()
        self.make_auto_sweep_gui ()
   
    def change_phase_scalarlimits(self, event=None):
        self.map.SetScalarRange((self.min_phcolor.get(),self.max_phcolor.get()))
	self.renwin.Render ()
 
    def do_normals (self, event=None):
        debug ("In IsoSurface::do_normals ()")
        Common.state.busy ()
        if self.normals_on.get ():
            self.norm.SetInput (self.cont_fil.GetOutput ())
            self.norm.SetFeatureAngle (self.angle_var.get ())
            self.map.SetInput (self.norm.GetOutput ())
        else:
            self.map.SetInput (self.cont_fil.GetOutput ())
        self.renwin.Render ()
        Common.state.idle ()

    def change_contour_limits (self, event=None):
        debug ("In IsoSurface::change_contour_limits ()")
        min_cnt = self.min_cnt.get ()
        max_cnt = self.max_cnt.get ()
        val = self.contour_var.get ()
        dr = self.data_range

        if max_cnt < val:
            msg = "Error: max. contour value less than current "\
                  "contour value."
            debug (msg)
            max_cnt = dr[1]
            self.max_cnt.set (max_cnt)

        if min_cnt > val:
            msg = "Error: min. contour value larger than current "\
                  "contour value."
            debug (msg)
            min_cnt = dr[0]
            self.min_cnt.set (min_cnt)

        resolution = (max_cnt - min_cnt)/1000
        if self.slider:
            self.slider.config (from_=min_cnt, to=max_cnt,
                                resolution=resolution)

    def change_contour_slider (self, event=None):
        debug ("In IsoSurface::change_contour_slider ()")
        Common.state.busy ()
        self.slider_pos = self.contour_var.get ()
        self.c_entry_var.set (self.slider_pos)
        self.cont_fil.SetValue (0, self.slider_pos)
        self.renwin.Render ()
        Common.state.idle ()

    def change_contour_entry (self, event=None):
        debug ("In IsoSurface::change_contour_entry ()")
        Common.state.busy ()
        self.slider_pos = self.c_entry_var.get ()
        self.contour_var.set (self.slider_pos)
        self.cont_fil.SetValue (0, self.slider_pos)
        self.renwin.Render ()
        Common.state.idle ()

    def do_sweep (self, event=None):
        debug ("In IsoSurface::do_sweep ()")
        if self.sweep_var.get ():
            val = int (1000*self.sweep_delay.get ())
            self.root.after (val, self.update_sweep)

    def update_sweep (self):
        debug ("In IsoSurface::update_sweep ()")
        if self.sweep_var.get ():
            min_cnt = self.min_cnt.get ()
            max_cnt = self.max_cnt.get ()
            d_pos = (max_cnt - min_cnt)/self.sweep_step.get ()
            pos = self.slider_pos + d_pos
            if (d_pos > 0) and (pos >= max_cnt):
                pos = min_cnt
            elif (d_pos < 0) and (pos <= min_cnt):
                pos = max_cnt
            self.contour_var.set (pos)
            self.change_contour_slider ()            
            val = int (1000*self.sweep_delay.get ())
            self.root.after (val, self.update_sweep)

    def close_gui (self, event=None):
        debug ("In IsoSurface::close_gui ()")
        Base.Objects.Module.close_gui (self, event)
        self.slider = None
