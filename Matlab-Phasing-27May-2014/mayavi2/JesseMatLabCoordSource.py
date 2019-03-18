from PPCoordSource import PPCoordSource
from scipy.io import loadmat
import numpy
import scipy.ndimage.measurements



class JesseMatLabCoordSource(PPCoordSource):

    def _open(self, fn):
        self._point_scalars_list.append(self.file_name)
        if not self.engine.dataarrays.has_key(self.file_name):
          f=loadmat(self.file_name)
	  array = f['array']
	  array = numpy.rot90(array,3)
	  
	  
	  nn=array.shape
	  cent=numpy.round((nn-numpy.array([1,1,1]))*.5)
	  
	  com=(numpy.round(scipy.ndimage.measurements.center_of_mass(numpy.abs(array))))
	  
	
	  
	  #temp=numpy.abs(array)
	  #com=(numpy.round( numpy.where(temp == temp.max())  ))
	    
	  print com
	  
	  #s0=numpy.int(cent[0]-com[0])
	  #s1=numpy.int(cent[1]-com[1])
	  #s2=numpy.int(cent[2]-com[2])
	  s0=0
	  s1=0
	  s2=0
	  
	  
	  array=numpy.roll(numpy.roll(numpy.roll(array,s0,axis=0),s1,axis=1),s2,axis=2)	
	  print 'centering data....'	  
	  
          self.engine.dataarrays[self.file_name]=array
          print "added new data to engine"
          del f
        else:
          print "can get data from engine"
        self.dataarray=self.engine.dataarrays[self.file_name]

