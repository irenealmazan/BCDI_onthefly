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



for file0 in os.listdir('.'): 
    if fnmatch.fnmatch(file0, '*AMP.rec'): 
        filename=file0 
	print filename 
        savename=filename	

	

	##
	scene1 = engine.new_scene()
	scene1.scene.off_screen_rendering = True

	jesse_mat_lab_coord_source1 = engine.open(filename, scene1)
	
	##
	from enthought.mayavi.modules.iso_surface import IsoSurface
	iso_surface1 = IsoSurface()
	engine.add_filter(iso_surface1, jesse_mat_lab_coord_source1)
	iso_surface1.contour.contours[0:1] = [0.5]
	
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
	    
	filename=dir+'/pics/AMP-'
	
	##
	module_manager1 = engine.scenes[1].children[0].children[0]
	module_manager1.scalar_lut_manager.show_scalar_bar = True
	module_manager1.scalar_lut_manager.show_legend = True
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.scalar_bar.number_of_labels = 3
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.number_of_labels = 3
	module_manager1.scalar_lut_manager.use_default_name = False
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.scalar_bar.title = ''
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.data_name = u''
	module_manager1.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
	module_manager1.scalar_lut_manager.label_text_property.italic = False
	module_manager1.scalar_lut_manager.label_text_property.font_family = 'arial'
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.data_range = array([ 0.4,  0.6 ])
	module_manager1.scalar_lut_manager.show_scalar_bar = False
	module_manager1.scalar_lut_manager.show_legend = False
	
	scene1.scene.x_plus_view()
	enthought.mayavi.mlab.savefig(filename+'+X.png', size=[1500,1500])
	enthought.mayavi.mlab.savefig(filename+'+X.ps', size=[1500,1500])
	scene1.scene.x_minus_view()
	enthought.mayavi.mlab.savefig(filename+'-X.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'-X.ps', size=[500,500])
	scene1.scene.y_plus_view()
	enthought.mayavi.mlab.savefig(filename+'+Y.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'+Y.ps', size=[500,500])
	scene1.scene.y_minus_view()
	enthought.mayavi.mlab.savefig(filename+'-Y.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'-Y.ps', size=[500,500])
	scene1.scene.z_plus_view()
	enthought.mayavi.mlab.savefig(filename+'+Z.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'+Z.ps', size=[500,500])
	scene1.scene.z_minus_view()
	enthought.mayavi.mlab.savefig(filename+'-Z.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'-Z.ps', size=[500,500])
		
	
	
	iso_surface1.contour.auto_contours = True
	iso_surface1.contour.minimum_contour = 0.4000000000000001
	iso_surface1.contour.maximum_contour = 0.6000000000000001
	iso_surface1.contour.number_of_contours = 3
	iso_surface1.actor.property.opacity = 0.29999999999999999

	scene1.scene.y_minus_view()
	scene1.scene.camera.elevation(15)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T+Y.png', size=[500,500])	
	enthought.mayavi.mlab.savefig(filename+'-T+Y.ps', size=[500,500])	

	scene1.scene.isometric_view()
	scene1.scene.camera.elevation(-25)
	scene1.scene.camera.azimuth(-5)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T-ISO.png', size=[500,500])	
	enthought.mayavi.mlab.savefig(filename+'-T-ISO.ps', size=[500,500])	

	module_manager1.scalar_lut_manager.show_scalar_bar = True
	
	#module_manager1 = engine.scenes[1].children[0].children[0]
	#module_manager1.scalar_lut_manager.scalar_bar_representation.minimum_size = array([1, 1])
	#module_manager1.scalar_lut_manager.scalar_bar_representation.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar_representation.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.scalar_bar_representation.maximum_size = array([100000, 100000])
	#module_manager1.scalar_lut_manager.scalar_bar.height = 0.80000000000000004
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.scalar_bar.width = 0.17000000000000001
	#module_manager1.scalar_lut_manager.show_scalar_bar = True
	#module_manager1.scalar_lut_manager.show_legend = True
	#module_manager1.scalar_lut_manager.use_default_name = False
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.scalar_bar.title = 'Test'
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.data_name = u''
	#module_manager1.scalar_lut_manager.use_default_range = False
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.data_range = array([ 0.2,  1. ])
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.data_range = array([ 0.2,  0.6])
	
	
	scene1.scene.y_minus_view()
	scene1.scene.camera.elevation(15)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T-SB+Y.png', size=[500,500])	
	enthought.mayavi.mlab.savefig(filename+'-T-SB+Y.ps', size=[500,500])

	scene1.scene.isometric_view()
	scene1.scene.camera.elevation(-25)
	scene1.scene.camera.azimuth(-5)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T-SB-ISO.png', size=[500,500])
	
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



for file0 in os.listdir('.'): 
    if fnmatch.fnmatch(file0, '*AMP-cjrf.rec'): 
        filename=file0 
	print filename 
        savename=filename	

	

	##
	scene1 = engine.new_scene()
	scene1.scene.off_screen_rendering = True

	jesse_mat_lab_coord_source1 = engine.open(filename, scene1)
	
	##
	from enthought.mayavi.modules.iso_surface import IsoSurface
	iso_surface1 = IsoSurface()
	engine.add_filter(iso_surface1, jesse_mat_lab_coord_source1)
	iso_surface1.contour.contours[0:1] = [0.5]
	
	##Orientation axes
	orientation_axes = OrientationAxes()
	engine.add_module(orientation_axes, obj=None)
	orientation_axes.text_property.shadow_offset = array([ 1, -1])
	orientation_axes.text_property.font_family = 'times'
	orientation_axes.text_property.shadow_offset = array([ 1, -1])
	orientation_axes.text_property.font_size = 15
	orientation_axes.axes.axis_labels = False
	##
	
	#if os.path.isdir(dir+'/pics/') == False:
	#    os.mkdir(dir+'/pics/')
	    
	filename=dir+'/pics/AMP-cjrf'
	
	##
	module_manager1 = engine.scenes[1].children[0].children[0]
	module_manager1.scalar_lut_manager.show_scalar_bar = True
	module_manager1.scalar_lut_manager.show_legend = True
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.scalar_bar.number_of_labels = 3
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.number_of_labels = 3
	module_manager1.scalar_lut_manager.use_default_name = False
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.scalar_bar.title = ''
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.data_name = u''
	module_manager1.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
	module_manager1.scalar_lut_manager.label_text_property.italic = False
	module_manager1.scalar_lut_manager.label_text_property.font_family = 'arial'
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager1.scalar_lut_manager.data_range = array([ 0.4,  0.6 ])
	module_manager1.scalar_lut_manager.show_scalar_bar = False
	module_manager1.scalar_lut_manager.show_legend = False
	
	scene1.scene.x_plus_view()
	enthought.mayavi.mlab.savefig(filename+'+X.png', size=[1500,1500])
	enthought.mayavi.mlab.savefig(filename+'+X.ps', size=[1500,1500])
	scene1.scene.x_minus_view()
	enthought.mayavi.mlab.savefig(filename+'-X.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'-X.ps', size=[500,500])
	scene1.scene.y_plus_view()
	enthought.mayavi.mlab.savefig(filename+'+Y.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'+Y.ps', size=[500,500])
	scene1.scene.y_minus_view()
	enthought.mayavi.mlab.savefig(filename+'-Y.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'-Y.ps', size=[500,500])
	scene1.scene.z_plus_view()
	enthought.mayavi.mlab.savefig(filename+'+Z.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'+Z.ps', size=[500,500])
	scene1.scene.z_minus_view()
	enthought.mayavi.mlab.savefig(filename+'-Z.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'-Z.ps', size=[500,500])
		
	
	
	iso_surface1.contour.auto_contours = True
	iso_surface1.contour.minimum_contour = 0.4000000000000001
	iso_surface1.contour.maximum_contour = 0.6000000000000001
	iso_surface1.contour.number_of_contours = 3
	iso_surface1.actor.property.opacity = 0.29999999999999999

	scene1.scene.y_minus_view()
	scene1.scene.camera.elevation(15)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T+Y.png', size=[500,500])	
	enthought.mayavi.mlab.savefig(filename+'-T+Y.ps', size=[500,500])	

	scene1.scene.isometric_view()
	scene1.scene.camera.elevation(-25)
	scene1.scene.camera.azimuth(-5)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T-ISO.png', size=[500,500])	
	enthought.mayavi.mlab.savefig(filename+'-T-ISO.ps', size=[500,500])	

	module_manager1.scalar_lut_manager.show_scalar_bar = True
	
	#module_manager1 = engine.scenes[1].children[0].children[0]
	#module_manager1.scalar_lut_manager.scalar_bar_representation.minimum_size = array([1, 1])
	#module_manager1.scalar_lut_manager.scalar_bar_representation.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar_representation.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.scalar_bar_representation.maximum_size = array([100000, 100000])
	#module_manager1.scalar_lut_manager.scalar_bar.height = 0.80000000000000004
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.scalar_bar.width = 0.17000000000000001
	#module_manager1.scalar_lut_manager.show_scalar_bar = True
	#module_manager1.scalar_lut_manager.show_legend = True
	#module_manager1.scalar_lut_manager.use_default_name = False
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.scalar_bar.title = 'Test'
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.data_name = u''
	#module_manager1.scalar_lut_manager.use_default_range = False
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.data_range = array([ 0.2,  1. ])
	#module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	#module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	#module_manager1.scalar_lut_manager.data_range = array([ 0.2,  0.6])
	
	
	scene1.scene.y_minus_view()
	scene1.scene.camera.elevation(15)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T-SB+Y.png', size=[500,500])	
	enthought.mayavi.mlab.savefig(filename+'-T-SB+Y.ps', size=[500,500])

	scene1.scene.isometric_view()
	scene1.scene.camera.elevation(-25)
	scene1.scene.camera.azimuth(-5)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T-SB-ISO.png', size=[500,500])
