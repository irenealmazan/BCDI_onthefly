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
print 'binx,biny:',binx, biny
det_px=temp['params']['det_px'] 
det_px=float(det_px[0][0][0][0]) 
det_py=det_px 
print 'det_px,py:',det_px, det_py
#arm = 0.55           #Detector distance from sample (units are m) 
arm=temp['params']['arm'] 
arm=float(arm[0][0][0][0]) 
print 'arm:', arm
lam=temp['params']['lam'] 
lam=float(lam[0][0][0][0]) 
print 'lam:', lam
#lam = 0.13933          #Wavelength of x-ray (units in nm) 
delta=temp['params']['delta'] 
delta=float(delta[0][0][0][0]) 
print 'delta:',delta
delta = delta * deg2rad    #Detector angles 
gam=temp['params']['gam'] 
gam=float(gam[0][0][0][0]) 
print 'gam:',gam
gam = gam * deg2rad    #Detector angles 
#gam =  13.6000 * deg2rad 
dth=temp['params']['dth'] 
dth=float(dth[0][0][0][0]) 
print 'dth:', dth
dth = dth * deg2rad    #Detector angles 
#dth = (1.0/50.0) * deg2rad      #step values for th translation 
dtilt=temp['params']['dtilt'] 
dtilt=float(dtilt[0][0][0][0]) 
dtilt = dtilt * deg2rad    #Detector angles 
#dtilt = 0.000 *deg2rad 
# load reconstructioed amp,phas and support for vtk output 

Qlabcenter1 = math.sin(delta)*math.cos(gam)*(2*3.14159265)/lam
Qlabcenter2 = math.sin(gam)*(2*3.14159265)/lam
Qlabcenter3 = (math.cos(delta)*math.cos(gam)-1.0)*(2*3.14159265)/lam
ki=(0.0,0.0,(2*3.14159265)/lam)
kf= ( math.sin(delta)*math.cos(gam)*(2*3.14159265)/lam, math.sin(gam)*(2*3.14159265)/lam, \
          math.cos(delta)*math.cos(gam)*(2*3.14159265)/lam )

vectorarr = vtk.vtkDoubleArray()
vectorarr.SetNumberOfComponents(3)
vectorarr.SetNumberOfTuples(3)
vectorarr.SetComponent(0,0,Qlabcenter1)
vectorarr.SetComponent(0,1,Qlabcenter2)
vectorarr.SetComponent(0,2,Qlabcenter3)
vectorarr.SetComponent(1,0,ki[0])
vectorarr.SetComponent(1,1,ki[1])
vectorarr.SetComponent(1,2,ki[2])
vectorarr.SetComponent(2,0,kf[0])
vectorarr.SetComponent(2,1,kf[1])
vectorarr.SetComponent(2,2,kf[2])

vectorpoints=vtk.vtkPoints()
vectorpoints.SetDataTypeToDouble()
vectorpoints.SetNumberOfPoints(3)
vectorpoints.SetPoint(0, (0.0,0.0,0.0) )
vectorpoints.SetPoint(1, (0.0,0.0,0.0) )
vectorpoints.SetPoint(2, (0.0,0.0,0.0) )

vectorgrid = vtk.vtkUnstructuredGrid()
vectorgrid.SetPoints(vectorpoints)
vectorgrid.GetPointData().SetVectors(vectorarr)

usgridwriter=vtk.vtkUnstructuredGridWriter()
usgridwriter.SetFileName('Qvector.vtk')
usgridwriter.SetFileTypeToASCII()
usgridwriter.SetInput(vectorgrid)
usgridwriter.Write()

###do the combining part, make a single vtk with both amplitude and phase

