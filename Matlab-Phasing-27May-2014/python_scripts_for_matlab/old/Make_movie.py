# Recorded script from Mayavi2
from numpy import array

import os
import fnmatch
from enthought.mayavi.modules.orientation_axes import OrientationAxes
import enthought.mayavi.mlab

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

	scene1 = engine.new_scene()
	scene1.scene.off_screen_rendering = True

	jesse_mat_lab_coord_source1 = engine.open(filename, scene1)
	from enthought.mayavi.modules.iso_surface import IsoSurface
	iso_surface1 = IsoSurface()
	engine.add_filter(iso_surface1, jesse_mat_lab_coord_source1)
	iso_surface1.contour.contours[0:1] = [0.5]
	
	filename='Movie-'
	os.mkdir(dir+'/Movie/')

	#scene1.scene.show_axes = True
	##Orientation axes
	orientation_axes = OrientationAxes()
	engine.add_module(orientation_axes, obj=None)
	orientation_axes.text_property.shadow_offset = array([ 1, -1])
	orientation_axes.text_property.font_family = 'times'
	orientation_axes.text_property.shadow_offset = array([ 1, -1])
	orientation_axes.text_property.font_size = 15
	orientation_axes.axes.axis_labels = False
	##
    
	

	int=5
	N=360/int

	for i in range(N):
	    
	    scene1.scene.x_plus_view()
	    scene1.scene.camera.azimuth(i*int)
	    scene1.scene.render()    
	    
	    #dir+'Movie/'+filename+str(i)+'+X.png'
	    name=dir+'/Movie/test'+str(i)+'.png'
	    
	    enthought.mayavi.mlab.savefig(name, size=[500,500])