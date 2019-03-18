import vtk
import sys

infile=sys.argv[1]
dimsin=sys.argv[2]

inbits=infile.split(".")
outputfile=inbits[0]+".bin"

sgr=vtk.vtkStructuredGridReader()
sgr.SetFileName(infile)
sgr.Update()


p_data = vtk.vtkStructuredPoints()
fil = vtk.vtkProbeFilter ()
fil.SetSource (sgr.GetOutput())
fil.SetInput(p_data)
fil.Update ()

pd = p_data
out = sgr.GetOutput()
b = out.GetBounds()
pd.SetOrigin(b[0], b[2], b[4])
l = [b[1] - b[0], b[3] - b[2], b[5] - b[4]]
tot_len = float(l[0] + l[1] + l[2])
npnt = pow(out.GetNumberOfPoints(), 1./3.) + 0.5
fac = 3.0*npnt/tot_len
dims=eval(dimsin)
pd.SetExtent(0, dims[0]-1, 0, dims[1] -1, 0, dims[2] -1)
pd.SetUpdateExtent(0, dims[0]-1, 0, dims[1] -1, 0, dims[2] -1)
pd.SetWholeExtent(0, dims[0]-1, 0, dims[1] -1, 0, dims[2] -1)
pd.SetDimensions(dims)
dims = [max(1, x-1) for x in dims]
l = [max(1e-3, x) for x in l]
sp = [l[0]/dims[0], l[1]/dims[1], l[2]/dims[2]]
pd.SetSpacing(sp)
pd.Update()

sp = pd.GetSpacing()
txt = "Point Spacing: [%.3g, %.3g, %.3g]"%(sp[0], sp[1], sp[2])
print txt
print p_data.GetDimensions()

#    def write_outputfile(self):
fil.Update()
pd=fil.GetOutput().GetPointData().GetScalars()
b=buffer(pd, 0, pd.GetDataTypeSize()*pd.GetSize())
f=open(outputfile, "wb")
f.write(b)
f.close()  
