# Recorded script from Mayavi2
from numpy import array

import os
import fnmatch
import sys

dir=os.getcwd()
os.chdir(dir)

save_dir='/Users/jesseclark/Documents/MATLAB/data_analysis/Au1111/J12-Al-100nm/T-analysis-6x10-CVL/vtks/movie_pics/combined/'

prefix='Ph-'
suffix='.jpg'
    
mod_n=int(sys.argv[1])
nthings=int(sys.argv[2])

#print mod_n
#print nthings

string='Ph*'

counter=1

for file0 in os.listdir('.'):
    if fnmatch.fnmatch(file0, string): 
    
    
        filename=file0 
	save_name=file0
    
	
    	rem=counter % nthings
	print rem
	if rem == mod_n:
	    cmd='cp '+dir+' '+filename+' '+save_dir
	    #print cmd
	    os.system(cmd)
	
	counter=counter+1
	
