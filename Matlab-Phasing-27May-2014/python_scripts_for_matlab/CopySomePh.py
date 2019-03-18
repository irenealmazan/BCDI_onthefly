# Recorded script from Mayavi2
from numpy import array

import os
import fnmatch

dir=os.getcwd()
os.chdir(dir)

base_dir=dir+'/'
save_dir=dir+'/combined/'

order=[1,7,2,3,6,4,5]

nthings=7

string='A-Ph*'

counter=1


for i in range(1,nthings+1):
	
	mod_n=int(order[i-1])

	if mod_n == nthings:
		mod_n=0

        print mod_n,i

        dir1=base_dir+'Ph'+str(i)+'/'
	os.chdir(dir1)
	print 'Ph'+str(i)

	counter = 1

	for file0 in os.listdir('.'):
	    #print file0
	    if fnmatch.fnmatch(file0, string): 
    
        	filename=file0 
		save_name=file0
    
   	 	rem=counter % nthings
	
		if rem == mod_n:
	    		cmd='cp '+dir1+filename+' '+save_dir
	    		print cmd
	    		os.system(cmd)
			
		counter=counter+1
		
	
