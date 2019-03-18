import vtk
#import mayavi as m
import sys
import os
import fnmatch

path=os.getcwd()
os.chdir(path)

for file0 in os.listdir('.'): 
    if fnmatch.fnmatch(file0, '*AMP.rec.vtk'): 
        filename=file0 
        basefilename=dir
        
        ampreader=vtk.vtkStructuredGridReader()
        ampreader.SetFileName(filename)
        ampreader.Update()
        
for file0 in os.listdir('.'): 
    if fnmatch.fnmatch(file0, '*PH.rec.vtk'):
        filename=file0 
        phreader=vtk.vtkStructuredGridReader()
        phreader.SetFileName(filename)
        phreader.Update()
        pharr=phreader.GetOutput().GetPointData().GetScalars()
        pharr.SetName("phases")
        
        data=ampreader.GetOutput()
        data.GetPointData().GetScalars().SetName("amps")
        data.GetPointData().AddArray( pharr )
        #data.GetPointData().SetActiveAttribute("phases", 0)
        data.Update()
        
        print file0
        filename=path+'/Amp-Phase.vtk'
        print filename
        gridwriter=vtk.vtkStructuredGridWriter()
        gridwriter.SetFileName(filename)
        gridwriter.SetFileTypeToBinary()
        gridwriter.SetInput(data)
        gridwriter.Write()

#v=m.mayavi()
#v.open_vtk_data( data )

#v.master.wait_window()
