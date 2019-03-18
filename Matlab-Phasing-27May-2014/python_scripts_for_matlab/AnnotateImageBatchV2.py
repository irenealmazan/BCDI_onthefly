from numpy import array

import os
import fnmatch

dir=os.getcwd()
os.chdir(dir)

base_dir=dir+'/'

prefix='Ph'

string_f='Ph*'

text=['030','030','107','155','155','201','235']

nthings=len(text)

for i in range(0,nthings):
    
    string=base_dir
    
    dir1=base_dir+prefix+str(i+1)+'/'
    os.chdir(dir1)

    

    for file0 in os.listdir('.'):
	
	if fnmatch.fnmatch(file0, string_f): 
    
	    filename=file0 
	    save_name='A-'+file0
	    
	    print filename
	    
	    text1='T='+text[i]+' c'
	    
	    cmd='convert '+filename+' -font ''Arial'' -pointsize 24 -gravity northwest -annotate 0 '+"'"+text1+"'"+' '+save_name
	    
	    os.system(cmd)