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
#in degrees
delta=32.5
gam=11.75
lam=.1377

##############################################
#convert to radians
delta=delta*math.pi/180
gam=gam*math.pi/180
########################
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
	
    
