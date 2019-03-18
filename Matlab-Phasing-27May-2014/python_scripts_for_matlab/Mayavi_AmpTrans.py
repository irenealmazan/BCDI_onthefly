# Recorded script from Mayavi2
from numpy import array

import os
import fnmatch
from enthought.mayavi.modules.orientation_axes import OrientationAxes
import enthought.mayavi.mlab
from enthought.mayavi.api import OffScreenEngine


try:
    engine = mayavi.engine
except NameError:
    from enthought.mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
# ------------------------------------------- 
dir=os.getcwd()

thresh=0.4
mag=1.
zoom=3.0
ph_max=3.14
ph_min=-3.14

rib=1

for file0 in os.listdir('.'): 
    if fnmatch.fnmatch(file0, 'Amp-Phase.vtk'): 
        filename=file0 
	print filename 
        savename=filename

	

	##
	scene = engine.new_scene()
	scene.scene.off_screen_rendering = True

	vtk_file_reader = engine.open(filename, scene)
	
	##
	
	
	##Orientation axes
	orientation_axes = OrientationAxes()
	engine.add_module(orientation_axes, obj=None)
	orientation_axes.text_property.shadow_offset = array([ 1, -1])
	orientation_axes.text_property.font_family = 'times'
	orientation_axes.text_property.shadow_offset = array([ 1, -1])
	orientation_axes.text_property.font_size = 15
	orientation_axes.axes.axis_labels = False
	##
	
	if os.path.isdir(dir+'/pics/') == False:
	    os.mkdir(dir+'/pics/')
	    
	filename=dir+'/pics/T-Iso-'
	
	##
	from enthought.mayavi.modules.iso_surface import IsoSurface
	iso_surface = IsoSurface()
	engine.add_filter(iso_surface, vtk_file_reader)
	iso_surface.actor.property.opacity = 0.29999999999999999
	iso_surface.contour.auto_contours = True
	iso_surface.contour.number_of_contours = 3
	iso_surface.contour.minimum_contour = thresh-.1
	iso_surface.contour.maximum_contour = thresh+.1
	module_manager = engine.scenes[1].children[0].children[0]
	module_manager.scalar_lut_manager.scalar_bar_representation.minimum_size = array([1, 1])
	module_manager.scalar_lut_manager.scalar_bar_representation.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar_representation.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.scalar_bar_representation.maximum_size = array([100000, 100000])
	module_manager.scalar_lut_manager.scalar_bar.height = 0.80000000000000004
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.scalar_bar.width = 0.17000000000000001
	module_manager.scalar_lut_manager.show_scalar_bar = True
	module_manager.scalar_lut_manager.show_legend = True
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.scalar_bar.number_of_labels = 3
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.number_of_labels = 3
	module_manager.scalar_lut_manager.use_default_range = False
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.data_range = array([ 0.3,  1. ])
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.data_range = array([ thresh-.1,thresh+.1])
	module_manager.scalar_lut_manager.use_default_name = False
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.scalar_bar.title = ''
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.data_name = u''
	module_manager.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
	module_manager.scalar_lut_manager.label_text_property.italic = False

	
	
	scene.scene.x_plus_view()
	scene.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'+X.png', size=[1500,1500],magnification='auto')
	enthought.mayavi.mlab.savefig(filename+'+X.ps', size=[1500,1500],magnification='auto')
	
	#enthought.mayavi.mlab.savefig(filename+'+X.rib', size=[1500,1500],magnification='auto')
	
	scene.scene.x_minus_view()
	scene.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'-X.png', size=[1500,1500],magnification='auto')
	enthought.mayavi.mlab.savefig(filename+'-X.ps', size=[1500,1500],magnification='auto')
	
	#enthought.mayavi.mlab.savefig(filename+'-X.rib', size=[1500,1500])
	
	scene.scene.y_plus_view()
	scene.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'+Y.png', size=[1500,1500],magnification='auto')
	enthought.mayavi.mlab.savefig(filename+'+Y.ps', size=[1500,1500],magnification='auto')
	scene.scene.y_minus_view()
	scene.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'-Y.png', size=[1500,1500],magnification='auto')
	enthought.mayavi.mlab.savefig(filename+'-Y.ps', size=[1500,1500],magnification='auto')
	scene.scene.z_plus_view()
	scene.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'+Z.png', size=[1500,1500],magnification='auto')
	enthought.mayavi.mlab.savefig(filename+'+Z.ps', size=[1500,1500],magnification='auto')
	scene.scene.z_minus_view()
	scene.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'-Z.png', size=[1500,1500],magnification='auto')
	enthought.mayavi.mlab.savefig(filename+'-Z.ps', size=[1500,1500],magnification='auto')
		

