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

thresh=0.25
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
	from enthought.mayavi.filters.contour import Contour
        contour = Contour()
        engine.add_filter(contour, obj=None)
        from enthought.mayavi.filters.set_active_attribute import SetActiveAttribute
        set_active_attribute = SetActiveAttribute()
        engine.add_filter(set_active_attribute, obj=None)
        from enthought.mayavi.modules.surface import Surface
        surface = Surface()
        engine.add_module(surface, obj=None)
	
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
	    
	filename=dir+'/pics/PhIso-fPi-'
	
	##
	module_manager = set_active_attribute.children[0]
        module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
        module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
        module_manager.scalar_lut_manager.data_range = array([-1.07110711,  1.11268517])
   	module_manager.scalar_lut_manager.default_data_range = array([-1.07110711,  1.11268517])
    	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.scalar_lut_manager.scalar_bar.title = 'Phase'
    	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.scalar_lut_manager.data_name = 'phases'
	module_manager.scalar_lut_manager.label_text_property.italic = False
	module_manager.scalar_lut_manager.label_text_property.font_family = 'arial'
    	module_manager.scalar_lut_manager.default_data_name = 'phases'
    	module_manager.vector_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.vector_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.vector_lut_manager.data_range = array([ 0.        ,  1.00000004])
    	module_manager.vector_lut_manager.default_data_range = array([ 0.        ,  1.00000004])
    	set_active_attribute.point_scalars_name = 'phases'
    	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.scalar_lut_manager.data_range = array([-1.07110711,  1.11268517])
    	module_manager.scalar_lut_manager.default_data_range = array([-1.07110711,  1.11268517])
	module_manager.scalar_lut_manager.data_range = array([ 0,  1 ])
    	module_manager.vector_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.vector_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.vector_lut_manager.data_range = array([ 0.,  1.])
    	module_manager.vector_lut_manager.default_data_range = array([ 0.,  1.])
    	module_manager.vector_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.vector_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.vector_lut_manager.scalar_bar.title = 'No vectors'
    	module_manager.vector_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.vector_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.vector_lut_manager.data_name = 'No vectors'
    	module_manager.vector_lut_manager.default_data_name = 'No vectors'
    	contour.filter.contours[0:1] = [thresh]
	module_manager.scalar_lut_manager.scalar_bar.title = ' '
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.data_name = u' '
	module_manager.scalar_lut_manager.data_range = array([ ph_min,  ph_max ])


	module_manager.scalar_lut_manager.show_scalar_bar = True
	module_manager.scalar_lut_manager.show_legend = True
	
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
		

ph_max=3.14/2
ph_min=-3.14/2

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
	from enthought.mayavi.filters.contour import Contour
        contour = Contour()
        engine.add_filter(contour, obj=None)
        from enthought.mayavi.filters.set_active_attribute import SetActiveAttribute
        set_active_attribute = SetActiveAttribute()
        engine.add_filter(set_active_attribute, obj=None)
        from enthought.mayavi.modules.surface import Surface
        surface = Surface()
        engine.add_module(surface, obj=None)
	
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
	    
	filename=dir+'/pics/PhIso-hPi-'
	
	##
	module_manager = set_active_attribute.children[0]
        module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
        module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
        module_manager.scalar_lut_manager.data_range = array([-1.07110711,  1.11268517])
   	module_manager.scalar_lut_manager.default_data_range = array([-1.07110711,  1.11268517])
    	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.scalar_lut_manager.scalar_bar.title = 'Phase'
    	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.scalar_lut_manager.data_name = 'phases'
	module_manager.scalar_lut_manager.label_text_property.italic = False
	module_manager.scalar_lut_manager.label_text_property.font_family = 'arial'
    	module_manager.scalar_lut_manager.default_data_name = 'phases'
    	module_manager.vector_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.vector_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.vector_lut_manager.data_range = array([ 0.        ,  1.00000004])
    	module_manager.vector_lut_manager.default_data_range = array([ 0.        ,  1.00000004])
    	set_active_attribute.point_scalars_name = 'phases'
    	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.scalar_lut_manager.data_range = array([-1.07110711,  1.11268517])
    	module_manager.scalar_lut_manager.default_data_range = array([-1.07110711,  1.11268517])
	module_manager.scalar_lut_manager.data_range = array([ 0,  1 ])
    	module_manager.vector_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.vector_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.vector_lut_manager.data_range = array([ 0.,  1.])
    	module_manager.vector_lut_manager.default_data_range = array([ 0.,  1.])
    	module_manager.vector_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.vector_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.vector_lut_manager.scalar_bar.title = 'No vectors'
    	module_manager.vector_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
    	module_manager.vector_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
    	module_manager.vector_lut_manager.data_name = 'No vectors'
    	module_manager.vector_lut_manager.default_data_name = 'No vectors'
    	contour.filter.contours[0:1] = [thresh]
	module_manager.scalar_lut_manager.scalar_bar.title = ' '
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.data_name = u' '
	module_manager.scalar_lut_manager.data_range = array([ ph_min,  ph_max ])


	module_manager.scalar_lut_manager.show_scalar_bar = True
	module_manager.scalar_lut_manager.show_legend = True
	
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
	

