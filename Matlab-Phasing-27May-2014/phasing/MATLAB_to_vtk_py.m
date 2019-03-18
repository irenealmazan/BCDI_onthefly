function MATLAB_to_vtk_py(dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% open a file for writing
fid = fopen([dir,'MATLAB_to_vtk.py'], 'w');

%% print a title, followed by a blank line
fprintf(fid,'import math \n');

fprintf(fid,'import math \n');
fprintf(fid,'import cstuff as c \n');
fprintf(fid,'import sys \n');
fprintf(fid,'import SpeFile as sf \n');
fprintf(fid,'import string as s \n');
fprintf(fid,'import vtk \n');
fprintf(fid,'import vtkImageImportFromArray as iia \n');
fprintf(fid,'import vtkStructuredGridMaker as vs \n');
fprintf(fid,'import vtkImageImportFromArray as iia \n');
fprintf(fid,'import fnmatch \n');
fprintf(fid,'import os \n');
fprintf(fid,'import numpy \n');
fprintf(fid,'from scipy.io import loadmat \n');

fprintf(fid,'################################################### \n');
fprintf(fid,'# load paramter file to get binning \n');

str=sprintf('for parf in os.listdir(''.''): \n');
fprintf(fid,str);
%fprintf(fid,'for parf in os.listdir('.'): \n');

str=sprintf('    if fnmatch.fnmatch(parf, ''*PARAMS.mat''): \n');
%fprintf(fid,'    if fnmatch.fnmatch(parf, '*PARAMS.mat'): \n');
fprintf(fid,str);

fprintf(fid,'        pfile=parf \n');
fprintf(fid,'	print pfile \n');

fprintf(fid,'temp =loadmat(pfile,struct_as_record=True)  \n');

fprintf(fid,'# Constants \n');
fprintf(fid,'pi=math.pi \n');
fprintf(fid,'deg2rad = pi/180.0000 \n');

str=sprintf('bin= temp[''params''][''binning''] \n');
fprintf(fid,str);
%fprintf(fid,'bin= temp['params']['binning'] \n');

fprintf(fid,'binx=float( bin[0][0][0][0]) \n');
fprintf(fid,'biny=float( bin[0][0][0][1]) \n');

str=sprintf('det_px=temp[''params''][''det_px''] \n');
fprintf(fid,str);
%fprintf(fid,'det_px=temp['params']['det_px'] \n');


fprintf(fid,'det_px=float(det_px[0][0][0][0]) \n');
fprintf(fid,'det_py=det_px \n');

fprintf(fid,'#arm = 0.55           #Detector distance from sample (units are m) \n');

str=sprintf('arm=temp[''params''][''arm''] \n');
%fprintf(fid,'arm=temp['params']['arm'] \n');
fprintf(fid,str);

fprintf(fid,'arm=float(arm[0][0][0][0]) \n');

str=sprintf('lam=temp[''params''][''lam''] \n');
fprintf(fid,str);
%fprintf(fid,'lam=temp['params']['lam'] \n');


fprintf(fid,'lam=float(lam[0][0][0][0]) \n');
fprintf(fid,'#lam = 0.13933          #Wavelength of x-ray (units in nm) \n');

str=sprintf('delta=temp[''params''][''delta''] \n');
fprintf(fid,str);
%fprintf(fid,'delta=temp['params']['delta'] \n');


fprintf(fid,'delta=float(delta[0][0][0][0]) \n');
fprintf(fid,'delta = delta * deg2rad    #Detector angles \n');

str=sprintf('gam=temp[''params''][''gam''] \n');
fprintf(fid,str);
%fprintf(fid,'gam=temp['params']['gam'] \n');


fprintf(fid,'gam=float(gam[0][0][0][0]) \n');
fprintf(fid,'gam = gam * deg2rad    #Detector angles \n');
fprintf(fid,'#gam =  13.6000 * deg2rad \n');

str=sprintf('dth=temp[''params''][''dth''] \n');
fprintf(fid,str);
%fprintf(fid,'dth=temp['params']['dth'] \n');


fprintf(fid,'dth=float(dth[0][0][0][0]) \n');
fprintf(fid,'dth = dth * deg2rad    #Detector angles \n');
fprintf(fid,'#dth = (1.0/50.0) * deg2rad      #step values for th translation \n');

str=sprintf('dtilt=temp[''params''][''dtilt''] \n');
fprintf(fid,str);
%fprintf(fid,'dtilt=temp['params']['dtilt'] \n');


fprintf(fid,'dtilt=float(dtilt[0][0][0][0]) \n');
fprintf(fid,'dtilt = dtilt * deg2rad    #Detector angles \n');

fprintf(fid,'#dtilt = 0.000 *deg2rad \n');

fprintf(fid,'# load reconstructioed amp,phas and support for vtk output \n');

str=sprintf('for file0 in os.listdir(''.''): \n');
fprintf(fid,str);
%fprintf(fid,'for file0 in os.listdir('.'): \n');

str=sprintf('    if fnmatch.fnmatch(file0, ''*.rec''): \n');
fprintf(fid,str);
%fprintf(fid,'    if fnmatch.fnmatch(file0, '*.rec'): \n');


fprintf(fid,'        filename=file0 \n');
fprintf(fid,'	print filename \n');

fprintf(fid,'	file=loadmat(filename) \n');

str=sprintf('	data=file[''array''] \n');
fprintf(fid,str);
%fprintf(fid,'	data=file['array'] \n');

fprintf(fid,'	data=numpy.rot90(data,3) \n');
fprintf(fid,'	size=data.shape \n');
	
fprintf(fid,'	#det_px=20.5e-6 \n');
fprintf(fid,'	#det_py=20.5e-6 \n');
	
fprintf(fid,'	print size \n');
fprintf(fid,'	# Put in array sizes and then what you want extracted \n');
fprintf(fid,'	size1=size[0] \n');
fprintf(fid,'	size2=size[1] \n');
fprintf(fid,'	size3=size[2] \n');
	
fprintf(fid,'	extractX = size1/2+10 		#xpositive-xnegative+10 \n');
fprintf(fid,'	extractY = size2/2+10 		#ypositive-ynegative+10 \n');
fprintf(fid,'	extractZ = size3/2+10		#zpositive-znegative+10 \n');
	
fprintf(fid,'	# Detector and scan variables \n');
fprintf(fid,'	pixelx=binx*(det_px)  #Total pixel size including binning (in meters) \n');
fprintf(fid,'	pixely=biny*(det_py) \n');
		
fprintf(fid,'	##################################################### \n');
	
fprintf(fid,'	diffData=c.Sp4Array() \n');
fprintf(fid,'	dDcrop=c.Sp4Array() \n');
fprintf(fid,'	supcrop=c.Sp4Array() \n');
fprintf(fid,'	support=c.Sp4Array() \n');
fprintf(fid,'	NewCoords=c.Sp4Array() \n');
	
fprintf(fid,'	#angles in RADIANS \n');
fprintf(fid,'	RotX=000.0*pi/180. \n');
fprintf(fid,'	RotY=000.0*pi/180. \n');
fprintf(fid,'	RotZ=000.0*pi/180. \n');
	
fprintf(fid,'	#################################################################### \n');
fprintf(fid,'	c.numpy_ToSp4(data,diffData) \n');
	
%fprintf(fid,'	#splitinputfilename=s.split(filename,'.') \n');
%fprintf(fid,'	#File=sf.SpeFile(filename) \n');
%fprintf(fid,'	#diffData=File.GetSp4Array() \n');
%fprintf(fid,'	#c.Sp4ArraySave(origarr,splitinputfilename[0]+'.sp4') \n');
fprintf(fid,'	#################################################################### \n');
	
fprintf(fid,'	#computed for use in the coordinate transformation \n');
fprintf(fid,'	dpy = pixely/arm #integer before * relates to binning i.e 2x2 is 2* \n');
fprintf(fid,'	dpx = pixelx/arm # "" \n');
	
fprintf(fid,'	if extractX>size1: \n');
fprintf(fid,'		extractX=size1 \n');
fprintf(fid,'	if extractY>size2: \n');
fprintf(fid,'		extractY=size2 \n');
fprintf(fid,'	if extractZ>size3: \n');
fprintf(fid,'		extractZ=size3 \n');
	
fprintf(fid,'	#Now lets produce vtk around just the support +-10 \n');
fprintf(fid,'	#diffData=seqdata.dist \n');
	
fprintf(fid,'	#Extract an array from the centre of the array which is defined by the sizes of Extractx,y,z inputted 	initially \n');
fprintf(fid,'	c.ExtractArray(diffData, dDcrop, size1/2-extractX/2, extractX, size2/2-extractY/2, extractY, size3/2-extractZ/2, extractZ) \n');
fprintf(fid,'	#c.ExtractArray(support, supcrop, size1/2-extractX/2, extractX, size2/2-extractY/2, extractY,size3/2-extractZ/2, extractZ) \n');
	
fprintf(fid,'	size1 = c.Sp4ArrayGetDim(dDcrop,0)[1] \n');
fprintf(fid,'	size2 = c.Sp4ArrayGetDim(dDcrop,1)[1] \n');
fprintf(fid,'	size3 = c.Sp4ArrayGetDim(dDcrop,2)[1] \n');
	
fprintf(fid,'	#define the dimensions in real space \n');
fprintf(fid,'	dx = 1.0/c.Sp4ArrayGetDim(diffData,0)[1] \n');
fprintf(fid,'	dy = 1.0/c.Sp4ArrayGetDim(diffData,1)[1] \n');
fprintf(fid,'	dz = 1.0/c.Sp4ArrayGetDim(diffData,2)[1] \n');
	
fprintf(fid,'	#Translate coordinates \n');
fprintf(fid,'	if dtilt ==0.0: \n');
fprintf(fid,'		 c.ThCoordTrans(NewCoords, size1, size2, size3, lam, delta, gam, dpx, dpy, dth,dx ,dy,dz,c.DIRECT) \n');
fprintf(fid,'		 print "theta rocking curve coordinate tranform" \n');
fprintf(fid,'	elif dth == 0.0: \n');
fprintf(fid,'		 c.TiltCoordTrans(NewCoords, size1, size2, size3, lam, delta, gam, dpx, dpy, dtilt,dx,dy,dz,c.DIRECT) \n');
fprintf(fid,'		 print "chi rocking curve coordinate transform" \n');
fprintf(fid,'	else: \n');
fprintf(fid,'		 print "dth/dtilt has been set incorrectly" \n');
	
	
fprintf(fid,'	#Ross added the Q-vector output 3/08 \n');
fprintf(fid,'	c.ShiftCoordOrigin(NewCoords, size1/2, size2/2, size3/2) \n');
	
fprintf(fid,'	s1=filename \n');%#		#seqdata.filenameroot+"dist%04d"%seqdata.iterationscompleted
	
%fprintf(fid,'	#s='%d_Support_%dPC-HIO_%dPO-ER'%(SequenceNumber,numHIOiter,numERiter) \n');
%fprintf(fid,'	#s=seqdata.filenameroot+"_Support" \n');
	
fprintf(fid,'	#vtk output \n');
fprintf(fid,'	dDampsGM=vs.vtkStructuredGridMaker(NewCoords, dDcrop, c.REAL)	#c.AMP \n');

str=sprintf('	dDampsGM.WriteStructuredGrid(s1+''.vtk'') \n');
fprintf(fid,str);
%fprintf(fid,'	dDampsGM.WriteStructuredGrid(s1+'.vtk') \n');
	
%fprintf(fid,'	#dDphaseGM=vs.vtkStructuredGridMaker(NewCoords, dDcrop, c.PHASE) \n');
%fprintf(fid,'	#dDphaseGM.WriteStructuredGrid(s1+'_phase.vtk') \n');
	
%fprintf(fid,'	#supampsGM=vs.vtkStructuredGridMaker(NewCoords, supcrop, c.AMP) \n');
%fprintf(fid,'	#if ShrinkWrap: \n');
%fprintf(fid,'	#  supampsGM.WriteStructuredGrid(s+'%d_amp.vtk'%seqdata.iterationscompleted) \n');
%fprintf(fid,'	#else: \n');
%fprintf(fid,'	#supampsGM.WriteStructuredGrid(s1+'_amp.vtk') \n');
	
fprintf(fid,'	#del dDampsGM \n');
fprintf(fid,'	#del dDphaseGM \n');
fprintf(fid,'	#del supampsGM \n');
	
fprintf(fid,'	c.Sp4ArrayDestroy(dDcrop) \n');
fprintf(fid,'	c.Sp4ArrayDestroy(supcrop) \n');
fprintf(fid,'	c.Sp4ArrayDestroy(NewCoords) \n');
	




%%
fclose(fid);

end

