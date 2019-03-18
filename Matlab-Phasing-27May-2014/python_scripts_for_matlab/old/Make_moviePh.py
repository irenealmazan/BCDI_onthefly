#sudo easy_install appscript

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
os.chdir(dir)

thresh=0.3			#threshold for isosurface
mag=1.3
zoom=3.0
ph_max=3.14/2
ph_min=-3.14/2

show_cbar= True      #True or False to show the legend

step_size=3			#rotation step size, degrees
images=int(round(360/step_size))	#total number of images

time=4
fps=25	#Frames per second, too little and mpeg can't handle it.  other formats can.
cjrf=0	#Will look for files that are of the conjugated and reflected reconstruction

transparent=0		#will mkae a movie with transparent isosurfaces

string='Amp-Phase.vtk'
    
for file0 in os.listdir('.'): 
    if fnmatch.fnmatch(file0, string): 
    
        filename=file0 
	print filename 
        savename=filename
	
	scene = engine.new_scene()
	scene.scene.off_screen_rendering = True
	vtk_file_reader = engine.open(filename, scene)
	
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

	if os.path.isdir(dir+'/movie_pics/') == False:
	    os.mkdir(dir+'/movie_pics/')

	if os.path.isdir(dir+'/movie_pics/Ph/') == False:
	    os.mkdir(dir+'/movie_pics/Ph/')

	    
	filename=dir+'/movie_pics/Ph/Ph-'
	this_dir=dir+'/movie_pics/Ph/'
	
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


	module_manager.scalar_lut_manager.show_scalar_bar = show_cbar
	module_manager.scalar_lut_manager.show_legend = show_cbar
	
	scene.scene.x_plus_view()
	scene.scene.camera.zoom(mag)
	
	angle=90
	for qq in range(1,images):
	    numb=str(qq)
	    
	    while len(numb) < len(str(2*images)):
		numb='0'+numb
	    
	    save_name=filename+numb+'.jpg'
	
	    print save_name
	    #enthought.mayavi.mlab.savefig(save_name+'+.png', size=[1500,1500])
	    scene.scene.save(save_name,size=[500,500])
	    scene.scene.camera.azimuth(-step_size)
	    scene.scene.render()

	
	for qq in range(1,images):
	    numb=str(qq+images-1)
	    
	    while len(numb) < len(str(2*images)):
		numb='0'+numb
	   
	    save_name=filename+numb+'.jpg'
	
	    print save_name
	    #enthought.mayavi.mlab.savefig(save_name+'+.png', size=[1500,1500])
	    scene.scene.save(save_name,size=[500,500])
	    scene.scene.camera.elevation(-step_size)
	    scene.scene.camera.orthogonalize_view_up()
	    scene.scene.render()

#########
print '<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
print numb
save_name='video.mpg'	    
if cjrf != 0:
    save_name='video-cjrf.mpg'
if transparent != 0:
    save_name='T-'+save_name
#########
#cmd_mpg='ffmpeg -f image2 -r '+str(fps)+' -vcodec mjpeg -qscale 1'+' -i '+this_dir+'%0'+str(len(numb))+'dcjrf.jpg  '+this_dir+save_name
cmd_mpg='ffmpeg -f image2 -r '+str(fps)+' -vcodec mjpeg -qscale 1'+' -i '+this_dir+'Ph-%0'+str(len(numb))+'d.jpg  '+this_dir+save_name
#######
#cmd_avi='ffmpeg -f image2 -r '+str(fps)+' -b 600k -i '+this_dir+'%0'+str(len(numb))+'d.jpg  '+this_dir+'video.avi'
#cmd_mpg='ffmpeg -f image2 -r '+str(fps)+' -b 1000k'+' -i '+this_dir+'%0'+str(len(numb))+'d.jpg  '+this_dir+'video.mpg'

#os.system(cmd_avi)
os.system(cmd_mpg)