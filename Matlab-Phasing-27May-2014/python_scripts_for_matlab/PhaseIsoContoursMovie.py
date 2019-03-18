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
os.chdir(dir)

thlow=-.75
thhigh=.5
mag=1.25

show_cbar=False

step_size=3			#rotation step size, degrees
images=int(round(360/step_size))	#total number of images
fps=25

xsize=500
ysize=500

opacity=.5

string='*Amp-Phase.vtk'

save_mname='/images-1/'

sstring='PhIso'

ref_file=dir+'/Run-40-Amp-Phase.vtk'

#--------------------------------------------

thing=1

for file0 in os.listdir('.'): 
    if fnmatch.fnmatch(file0, string): 

        filename=file0
        
        if os.path.isdir(dir+save_mname) == False:
	    os.mkdir(dir+save_mname)
        
        if os.path.isdir(dir+save_mname+sstring+str(thing)+'/') == False:
	    os.mkdir(dir+save_mname+sstring+str(thing)+'/')

	this_dir=dir+save_mname+sstring+str(thing)+'/'    
	files_string=this_dir+sstring+'-'
        thing=thing+1


        scene = engine.scenes[0]
        vtk_file_reader = engine.open(ref_file,scene)
        from enthought.mayavi.modules.iso_surface import IsoSurface
        iso_surface = IsoSurface()
        engine.add_filter(iso_surface, vtk_file_reader)
        iso_surface.contour.contours[0:1] = [0.050000000000000003]
        iso_surface.actor.property.opacity = 0.10000000000000001
        vtk_file_reader1 = engine.open(filename)
        vtk_file_reader1.point_scalars_name = 'phases'
        
        iso_surface1 = IsoSurface()
        engine.add_filter(iso_surface1, vtk_file_reader1)
        iso_surface1.contour.auto_contours = True
        iso_surface1.contour.number_of_contours = 2
        iso_surface1.contour.minimum_contour = thlow
        iso_surface1.contour.maximum_contour = thhigh
        module_manager1 = vtk_file_reader1.children[0]
        module_manager1.scalar_lut_manager.use_default_range = False
        module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
        module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
        module_manager1.scalar_lut_manager.data_range = array([-0.5     ,  2.610557])
        module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
        module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
        module_manager1.scalar_lut_manager.data_range = array([thlow,  thhigh])
        
	iso_surface1.actor.property.opacity = opacity

	##Orientation axes
	orientation_axes = OrientationAxes()
	engine.add_module(orientation_axes, obj=None)
	orientation_axes.text_property.shadow_offset = array([ 1, -1])
	orientation_axes.text_property.font_family = 'times'
	orientation_axes.text_property.shadow_offset = array([ 1, -1])
	orientation_axes.text_property.font_size = 15
	orientation_axes.axes.axis_labels = False
	
	scene.scene.y_minus_view()
	scene.scene.camera.zoom(mag)
	
	#color bar
	module_manager1.scalar_lut_manager.scalar_bar_representation.minimum_size = array([1, 1])
	module_manager1.scalar_lut_manager.scalar_bar_representation.position2 = array([ 0.17,  0.8 ])
	module_manager1.scalar_lut_manager.scalar_bar_representation.position = array([ 0.96,  0.1 ])
	module_manager1.scalar_lut_manager.scalar_bar_representation.maximum_size = array([100000, 100000])
	module_manager1.scalar_lut_manager.scalar_bar.height = 0.80000000000000004
	module_manager1.scalar_lut_manager.scalar_bar.position2 = array([ 0.4,  0.9 ])
	module_manager1.scalar_lut_manager.scalar_bar.position = array([ 0.9,  0.1 ])
	module_manager1.scalar_lut_manager.scalar_bar.width = 0.1000000000000001
	module_manager1.scalar_lut_manager.show_scalar_bar = show_cbar
	module_manager1.scalar_lut_manager.show_legend = show_cbar
	
        
	
	angle=90
	for qq in range(1,images):
	    numb=str(qq)
	    
	    while len(numb) < len(str(2*images)):
		numb='0'+numb
	    
	    save_name=files_string+numb+'.jpg'
	
	    print save_name
	    #enthought.mayavi.mlab.savefig(save_name+'+.png', size=[1500,1500])
	    scene.scene.save(save_name,size=[xsize,ysize])
	    scene.scene.camera.azimuth(-step_size)
	    scene.scene.render()
        
	for qq in range(1,images):
	    numb=str(qq+images-1)
	    
	    while len(numb) < len(str(2*images)):
		numb='0'+numb
	   
	    save_name=files_string+numb+'.jpg'
	
	    print save_name
	    #enthought.mayavi.mlab.savefig(save_name+'+.png', size=[1500,1500])
	    scene.scene.save(save_name,size=[xsize,ysize])
	    scene.scene.camera.elevation(-step_size)
	    scene.scene.camera.orthogonalize_view_up()
	    scene.scene.render()
	    
	save_name=sstring+'.mpg'	
	cmd_mpg='ffmpeg -f image2 -r '+str(fps)+' -vcodec mjpeg -qscale 1'+' -i '+this_dir+sstring+'-%0'+str(len(numb))+'d.jpg  '+this_dir+save_name
	os.system(cmd_mpg)