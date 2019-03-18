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
mag=1.3

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
	iso_surface1.contour.contours[0:1] = [thresh]
	
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
	module_manager1.scalar_lut_manager.data_range = array([ thresh-0.1,  thresh+0.1 ])
	module_manager1.scalar_lut_manager.show_scalar_bar = False
	module_manager1.scalar_lut_manager.show_legend = False
	
	scene1.scene.x_plus_view()
	scene1.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'+X.png', size=[1500,1500])
	enthought.mayavi.mlab.savefig(filename+'+X.ps', size=[1500,1500])
	scene1.scene.x_minus_view()
	scene1.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'-X.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'-X.ps', size=[500,500])
	scene1.scene.y_plus_view()
	scene1.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'+Y.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'+Y.ps', size=[500,500])
	scene1.scene.y_minus_view()
	scene1.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'-Y.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'-Y.ps', size=[500,500])
	scene1.scene.z_plus_view()
	scene1.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'+Z.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'+Z.ps', size=[500,500])
	scene1.scene.z_minus_view()
	scene1.scene.camera.zoom(mag)
	enthought.mayavi.mlab.savefig(filename+'-Z.png', size=[500,500])
	enthought.mayavi.mlab.savefig(filename+'-Z.ps', size=[500,500])
		
	#scalar cut plane
	module_manager1.scalar_lut_manager.data_range = array([ 0,  1 ])
	scene1.scene.x_plus_view()
	iso_surface1 = engine.scenes[1].children[0].children[0].children[0]
	iso_surface1.actor.property.opacity = 0.20000000000000001
	from enthought.mayavi.modules.scalar_cut_plane import ScalarCutPlane
	scalar_cut_plane = ScalarCutPlane()
	engine.add_module(scalar_cut_plane, obj=None)

	origin=scalar_cut_plane.implicit_plane.plane.origin
	print 'Origin -',origin

	scalar_cut_plane.implicit_plane.widget.origin = origin
	scalar_cut_plane.warp_scalar.filter.normal = array([ 1.,  0.,  0.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 1.,  0.,  0.])
	scalar_cut_plane.implicit_plane.widget.enabled = False
	scene1.scene.reset_zoom()
	scene1.scene.camera.zoom(3.5)
	scene1.scene.render()
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
	module_manager1.scalar_lut_manager.label_text_property.font_family = 'arial'
	module_manager1.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
	module_manager1.scalar_lut_manager.label_text_property.color = (1.0, 1.0, 1.0)
	
	enthought.mayavi.mlab.savefig(filename+'+SCP+X.png', size=[500,500])
	scene1.scene.y_plus_view()
	scalar_cut_plane.implicit_plane.widget.origin = origin
	scalar_cut_plane.warp_scalar.filter.normal = array([ 1.,  0.,  0.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 1.,  0.,  0.])
	scalar_cut_plane.implicit_plane.widget.origin = origin
	scalar_cut_plane.implicit_plane.plane.origin = origin
	scalar_cut_plane.implicit_plane.plane.global_warning_display = 0
	scalar_cut_plane.implicit_plane.plane.normal = array([ 0.,  1.,  0.])
	scalar_cut_plane.warp_scalar.filter.global_warning_display = 0
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  1.,  0.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  1.,  0.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 0.,  1.,  0.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  0.,  0.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  0.,  0.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 0.,  0.,  0.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  1.,  0.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  1.,  0.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 0.,  1.,  0.])
	scalar_cut_plane.implicit_plane.widget.normal_to_y_axis = True
	scene1.scene.reset_zoom()
	scene1.scene.camera.zoom(3.5)
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'+SCP+Y.png', size=[1500,1500])
	scene1.scene.z_plus_view()
	scalar_cut_plane.implicit_plane.widget.origin = origin
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  1.,  0.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 0.,  1.,  0.])
	scalar_cut_plane.implicit_plane.widget.normal_to_y_axis = 0
	scalar_cut_plane.implicit_plane.plane.origin = origin
	scalar_cut_plane.implicit_plane.plane.normal = array([ 0.,  0.,  1.])
	scalar_cut_plane.implicit_plane.widget.origin = origin
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  0.,  1.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  0.,  1.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 0.,  0.,  1.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  0.,  0.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  0.,  0.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 0.,  0.,  0.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  0.,  1.])
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  0.,  1.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 0.,  0.,  1.])
	scalar_cut_plane.implicit_plane.widget.normal_to_z_axis = True
	scene1.scene.reset_zoom()
	scene1.scene.camera.zoom(3.5)
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'+SCP+Z.png', size=[1500,1500])
	scalar_cut_plane.implicit_plane.widget.origin = origin
	scalar_cut_plane.warp_scalar.filter.normal = array([ 0.,  0.,  1.])
	scalar_cut_plane.implicit_plane.widget.normal = array([ 0.,  0.,  1.])
	scalar_cut_plane.implicit_plane.widget.interactor = None
	module_manager1 = engine.scenes[1].children[0].children[0]
	module_manager1.children[2:3] = []
	iso_surface1.actor.property.opacity = 1.0
	scene1.scene.x_plus_view()
	scene1.scene.reset_zoom()
	scene1.scene.camera.zoom(3.5)
	scene1.scene.render()
	module_manager1.scalar_lut_manager.show_scalar_bar = False
	module_manager1.scalar_lut_manager.show_legend = False
	module_manager1.scalar_lut_manager.data_range = array([ thresh-0.1,  thresh+0.1 ])
	module_manager1.scalar_lut_manager.label_text_property.color = (0, 0, 0)
	##Orientation axes
	#orientation_axes = OrientationAxes()
	#engine.add_module(orientation_axes, obj=None)
	#orientation_axes.text_property.shadow_offset = array([ 1, -1])
	#orientation_axes.text_property.font_family = 'times'
	#orientation_axes.text_property.shadow_offset = array([ 1, -1])
	#orientation_axes.text_property.font_size = 15
	#orientation_axes.axes.axis_labels = False
	##
	#Transparent isosurfaces
	iso_surface1.contour.auto_contours = True
	iso_surface1.contour.minimum_contour = thresh-0.1
	iso_surface1.contour.maximum_contour = thresh+0.1
	iso_surface1.contour.number_of_contours = 3
	iso_surface1.actor.property.opacity = 0.29999999999999999

	scene1.scene.y_minus_view()
	#scene1.scene.camera.zoom(3.5)
	#scene1.scene.render()
	scene1.scene.camera.elevation(15)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.camera.zoom(mag)
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T+Y.png', size=[500,500])	
	enthought.mayavi.mlab.savefig(filename+'-T+Y.ps', size=[500,500])	

	scene1.scene.isometric_view()
	scene1.scene.camera.elevation(-25)
	scene1.scene.camera.azimuth(-5)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.camera.zoom(mag)
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T-ISO.png', size=[500,500])	
	enthought.mayavi.mlab.savefig(filename+'-T-ISO.ps', size=[500,500])	

	module_manager1.scalar_lut_manager.show_scalar_bar = True
	
	
	scene1.scene.y_minus_view()
	scene1.scene.camera.elevation(15)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.camera.zoom(mag)
	scene1.scene.render()
	enthought.mayavi.mlab.savefig(filename+'-T-SB+Y.png', size=[500,500])	
	enthought.mayavi.mlab.savefig(filename+'-T-SB+Y.ps', size=[500,500])

	scene1.scene.isometric_view()
	scene1.scene.camera.elevation(-25)
	scene1.scene.camera.azimuth(-5)
	scene1.scene.camera.orthogonalize_view_up()
	scene1.scene.camera.zoom(mag)
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
	iso_surface1.contour.contours[0:1] = [thresh]
	
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
	module_manager1.scalar_lut_manager.data_range = array([ thresh-0.1,  thresh+0.1 ])
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
	iso_surface1.contour.minimum_contour = thresh-0.1
	iso_surface1.contour.maximum_contour = thresh+0.1
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
