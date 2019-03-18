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

thresh=0.25			#threshold for isosurface
mag=1.3
zoom=3.0
ph_max=1.5
ph_min=-1.5

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
	
	if os.path.isdir(dir+'/movie_pics/') == False:
	    os.mkdir(dir+'/movie_pics/')
	    
	filename=dir+'/movie_pics/'
	this_dir=filename
	
	from enthought.mayavi.modules.iso_surface import IsoSurface
	iso_surface = IsoSurface()
	engine.add_filter(iso_surface, obj=None)
	iso_surface.contour.contours[0:1] = [thresh]
	scene = engine.scenes[0]
	scene.scene.render()
	scene.scene.x_plus_view()
	scene.scene.camera.zoom(mag)
	#scene.scene.render()
	## change the color bar range
	module_manager = engine.scenes[0].children[0].children[0]
	module_manager.scalar_lut_manager.show_scalar_bar = True
	module_manager.scalar_lut_manager.show_legend = True
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.scalar_bar.number_of_labels = 3
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.number_of_labels = 3
	module_manager.scalar_lut_manager.use_default_name = False
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.scalar_bar.title = ''
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.data_name = u''
	module_manager.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
	module_manager.scalar_lut_manager.label_text_property.italic = False
	module_manager.scalar_lut_manager.label_text_property.font_family = 'arial'
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.17,  0.8 ])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.82,  0.1 ])
	module_manager.scalar_lut_manager.data_range = array([ thresh-0.1,  thresh+0.1 ])
	module_manager.scalar_lut_manager.show_scalar_bar = False
	module_manager.scalar_lut_manager.show_legend = False
	
	if transparent != 0:
	    iso_surface.contour.auto_contours = True
	    iso_surface.contour.minimum_contour = thresh-0.1
	    iso_surface.contour.maximum_contour = thresh+0.1
	    iso_surface.contour.number_of_contours = 3
	    iso_surface.actor.property.opacity = 0.29999999999999999
	
	angle=0
	for qq in range(1,images):
	    numb=str(qq)
	    
	    while len(numb) < len(str(images)):
		numb='0'+numb
	    #print numb,len(str(images))
	    
	    #if cjrf != 0:
		#save_name=filename+numb+'cjrf.jpg'
	    #else:
	    save_name=filename+numb+'.jpg'
	
	    print save_name
	    #enthought.mayavi.mlab.savefig(save_name+'+.png', size=[1500,1500])
	    scene.scene.save(save_name,size=[500,500])
	    scene.scene.camera.azimuth(-step_size)
	    scene.scene.render()
#########
save_name='video.mpg'	    
if cjrf != 0:
    save_name='video-cjrf.mpg'
if transparent != 0:
    save_name='T-'+save_name
#########
#cmd_mpg='ffmpeg -f image2 -r '+str(fps)+' -vcodec mjpeg -qscale 1'+' -i '+this_dir+'%0'+str(len(numb))+'dcjrf.jpg  '+this_dir+save_name
cmd_mpg='ffmpeg -f image2 -r '+str(fps)+' -vcodec mjpeg -qscale 1'+' -i '+this_dir+'%0'+str(len(numb))+'d.jpg  '+this_dir+save_name
#######
#cmd_avi='ffmpeg -f image2 -r '+str(fps)+' -b 600k -i '+this_dir+'%0'+str(len(numb))+'d.jpg  '+this_dir+'video.avi'
#cmd_mpg='ffmpeg -f image2 -r '+str(fps)+' -b 1000k'+' -i '+this_dir+'%0'+str(len(numb))+'d.jpg  '+this_dir+'video.mpg'

#os.system(cmd_avi)
os.system(cmd_mpg)