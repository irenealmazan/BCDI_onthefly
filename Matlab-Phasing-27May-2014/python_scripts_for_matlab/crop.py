from numpy import array

import os
import fnmatch

dir=os.getcwd()
os.chdir(dir)

base_dir=dir+'/'

string_f='*I-Ph-xy.png'

for file0 in os.listdir('.'):
	
	if fnmatch.fnmatch(file0, string_f): 
    
	    filename=file0 
	    save_name='A-'+file0
	    
	    print filename
	    
	    #cmd='convert '+filename+' -font ''Arial'' -pointsize 24 -gravity northwest -annotate 0 '+"'"+text1+"'"+' '+save_name
	    cmd='convert '+filename+' -crop 1582x826+348+488 '+save_name

	    os.system(cmd)