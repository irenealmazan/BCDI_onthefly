import fileinput
import sys

def dirEntries(dir_name, subdir, *args):
    '''Return a list of file names found in directory 'dir_name'
    If 'subdir' is True, recursively access subdirectories under 'dir_name'.
    Additional arguments, if any, are file extensions to match filenames. Matched
        file names are added to the list.
    If there are no additional arguments, all files found in the directory are
        added to the list.
    Example usage: fileList = dirEntries(r'H:\TEMP', False, 'txt', 'py')
        Only files with 'txt' and 'py' extensions will be added to the list.
    Example usage: fileList = dirEntries(r'H:\TEMP', True)
        All files and all the files in subdirectories under H:\TEMP will be added
        to the list.
    '''
    fileList = []
    for file in os.listdir(dir_name):
        dirfile = os.path.join(dir_name, file)
        if os.path.isfile(dirfile):
            if not args:
                fileList.append(dirfile)
            else:
                if os.path.splitext(dirfile)[1][1:] in args:
                    fileList.append(dirfile)
        # recursively access file names in subdirectories
        elif os.path.isdir(dirfile) and subdir:
            print "Accessing directory:", dirfile
            fileList.extend(dirEntries(dirfile, subdir, *args))
    return fileList


def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
                line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

import os
import re
import string

dir=os.getcwd()
filelist=dirEntries(dir,True,'py')

search_file='MATLAB_to_vtk_Qv.py'
#search_file='test_batch.py'

length=len(search_file)

for file0 in filelist:
	
        file_len=len(file0)
        
        thing=file0[file_len-length:file_len]
        
        check= (thing == search_file) 
        
        if string.find(file0, search_file) != -1:
            print file0
            
            cdir=file0[0:len(file0)-length]
            print cdir
            
            
            
            
            s2="os.listdir('"+cdir+"'):"
            #
            s3="os.listdir('.'):"
            #
            #
            replaceAll(file0,s3,s2)
            #
            s4="file=loadmat(filename)" 
            s5="file=loadmat("+"'"+cdir+"'+"+"filename)"
            #
            replaceAll(file0,s4,s5)
            #
            #
            #
            s6="s1=filename"
            s7="s1="+"'"+cdir+"'"+"+"+"filename"    
            replaceAll(file0,s6,s7)
            
	    s8="usgridwriter.SetFileName('Qvector.vtk')"
	    s9="usgridwriter.SetFileName('"+cdir+"Qvector.vtk')"
	    replaceAll(file0,s8,s9)

	    s10="path=os.getcwd()"
	    s11="path='"+cdir+"'"
	    replaceAll(file0,s10,s11)

	    
            os.system('python '+file0)
