# Recorded script from Mayavi2
from numpy import array

import os
import fnmatch

dir=os.getcwd()
os.chdir(dir)

prefix='Ph-'
suffix='.jpg'
    
text='t=+010 ps'

string='Ph*'

for file0 in os.listdir('.'): 
    if fnmatch.fnmatch(file0, string): 
    
        filename=file0 
	save_name='A-'+file0
	
	cmd='convert '+filename+' -font ''Arial'' -pointsize 24 -gravity northwest -annotate 0 '+"'"+text+"'"+' '+save_name
	
	os.system(cmd)
	
	
