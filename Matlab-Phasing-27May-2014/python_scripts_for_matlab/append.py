from numpy import array

import os
import fnmatch

dir=os.getcwd()
os.chdir(dir)

base_dir=dir+'/'

string_f='A-*Ph-zy*.png'

qq=0

savename='appended.png'

for file0 in os.listdir('.'):
	
	if fnmatch.fnmatch(file0, string_f): 
    
	    filename=file0 
	    save_name='A-'+file0
	    qq=qq+1	    
	    print filename
	    name2=filename
	    
           

	    if qq > 1:
            	cmd='convert '+name1+' '+name2+' +append '+savename

	    	os.system(cmd)
	    	name1=savename
 

	    if qq == 1:
		name1=filename

 	    