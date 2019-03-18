import math 
import math 
import cstuff as c 
import sys 
import SpeFile as sf 
import string as s 
import vtk 
import vtkImageImportFromArray as iia 
import vtkStructuredGridMaker as vs 
import vtkImageImportFromArray as iia 
import fnmatch 
import os 
import numpy 
from scipy.io import loadmat 
################################################### 
# load paramter file to get binning 
for parf in os.listdir('.'): 
    if fnmatch.fnmatch(parf, '*PARAMS.mat'): 
        pfile=parf 
	print pfile 
temp =loadmat(pfile,struct_as_record=True)  
# Constants 
pi=math.pi 
deg2rad = pi/180.0000 
bin= temp['params']['binning'] 
binx=float( bin[0][0][0][0]) 
biny=float( bin[0][0][0][1]) 
det_px=temp['params']['det_px'] 
det_px=float(det_px[0][0][0][0]) 
det_py=det_px 
#arm = 0.55           #Detector distance from sample (units are m) 
arm=temp['params']['arm'] 
arm=float(arm[0][0][0][0]) 
lam=temp['params']['lam'] 
lam=float(lam[0][0][0][0]) 
#lam = 0.13933          #Wavelength of x-ray (units in nm) 
delta=temp['params']['delta'] 
delta=float(delta[0][0][0][0]) 
delta = delta * deg2rad    #Detector angles 
gam=temp['params']['gam'] 
gam=float(gam[0][0][0][0]) 
gam = gam * deg2rad    #Detector angles 
#gam =  13.6000 * deg2rad 
dth=temp['params']['dth'] 
dth=float(dth[0][0][0][0]) 
dth = dth * deg2rad    #Detector angles 
#dth = (1.0/50.0) * deg2rad      #step values for th translation 
dtilt=temp['params']['dtilt'] 
dtilt=float(dtilt[0][0][0][0]) 
dtilt = dtilt * deg2rad    #Detector angles 
#dtilt = 0.000 *deg2rad 
# load reconstructioed amp,phas and support for vtk output 
for file0 in os.listdir('.'): 
    if fnmatch.fnmatch(file0, '*COH.rec'): 
        filename=file0 
	print filename 
	file=loadmat(filename) 
	data=file['array'] 
	data=numpy.rot90(data,3) 
	data=abs(data)
	size=data.shape 
	#det_px=20.5e-6 
	#det_py=20.5e-6 
	print size 
	# Put in array sizes and then what you want extracted 
	size1=size[0] 
	size2=size[1] 
	size3=size[2] 
	extractX = size1#/2+10 		#xpositive-xnegative+10 
	extractY = size2#/2+10 		#ypositive-ynegative+10 
	extractZ = size3#/2+10		#zpositive-znegative+10 
	# Detector and scan variables 
	pixelx=binx*(det_px)  #Total pixel size including binning (in meters) 
	pixely=biny*(det_py) 
	##################################################### 
	diffData=c.Sp4Array() 
	dDcrop=c.Sp4Array() 
	supcrop=c.Sp4Array() 
	support=c.Sp4Array() 
	NewCoords=c.Sp4Array() 
	#angles in RADIANS 
	RotX=000.0*pi/180. 
	RotY=000.0*pi/180. 
	RotZ=000.0*pi/180. 
	#################################################################### 
	c.numpy_ToSp4(data,diffData) 
	#################################################################### 
	#computed for use in the coordinate transformation 
	dpy = pixely/arm #integer before * relates to binning i.e 2x2 is 2* 
	dpx = pixelx/arm # "" 
	if extractX>size1: 
		extractX=size1 
	if extractY>size2: 
		extractY=size2 
	if extractZ>size3: 
		extractZ=size3 
	#Now lets produce vtk around just the support +-10 
	#diffData=seqdata.dist 
	#Extract an array from the centre of the array which is defined by the sizes of Extractx,y,z inputted 	initially 
	c.ExtractArray(diffData, dDcrop, size1/2-extractX/2, extractX, size2/2-extractY/2, extractY, size3/2-extractZ/2, extractZ) 
	#c.ExtractArray(support, supcrop, size1/2-extractX/2, extractX, size2/2-extractY/2, extractY,size3/2-extractZ/2, extractZ) 
	size1 = c.Sp4ArrayGetDim(dDcrop,0)[1] 
	size2 = c.Sp4ArrayGetDim(dDcrop,1)[1] 
	size3 = c.Sp4ArrayGetDim(dDcrop,2)[1] 
	#define the dimensions in real space 
	dx = 1.0/c.Sp4ArrayGetDim(diffData,0)[1] 
	dy = 1.0/c.Sp4ArrayGetDim(diffData,1)[1] 
	dz = 1.0/c.Sp4ArrayGetDim(diffData,2)[1] 
	#Translate coordinates 
	if dtilt ==0.0: 
		 c.ThCoordTrans(NewCoords, size1, size2, size3, lam, delta, gam, dpx, dpy, dth,dx ,dy,dz,c.DIRECT) 
		 print "theta rocking curve coordinate tranform" 
	elif dth == 0.0: 
		 c.TiltCoordTrans(NewCoords, size1, size2, size3, lam, delta, gam, dpx, dpy, dtilt,dx,dy,dz,c.DIRECT) 
		 print "chi rocking curve coordinate transform" 
	else: 
		 print "dth/dtilt has been set incorrectly" 
	#Ross added the Q-vector output 3/08 
	c.ShiftCoordOrigin(NewCoords, size1/2, size2/2, size3/2) 
	s1=filename 
	#vtk output 
	dDampsGM=vs.vtkStructuredGridMaker(NewCoords, dDcrop, c.AMP)	#c.AMP 
	dDampsGM.WriteStructuredGrid(s1+'.vtk') 
	#del dDampsGM 
	#del dDphaseGM 
	#del supampsGM 
	c.Sp4ArrayDestroy(dDcrop) 
	c.Sp4ArrayDestroy(supcrop) 
	c.Sp4ArrayDestroy(NewCoords) 
