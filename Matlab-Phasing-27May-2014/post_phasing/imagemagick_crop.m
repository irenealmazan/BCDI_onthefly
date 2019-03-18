function [ output_args ] = imagemagick_crop( input_args )
%jclark
%use imagemagick to crop images

%recursive string
string_f='*Ph-zy*.png';

%get current dir
name_of_this_file='imagemagick_crop';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
fdir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

%search for the file names
fnames=rdir([fdir,'/**/*/',string_f]);


        filename=file0 
	    save_name='A-'+file0
	    
	    print filename
	    
	    #cmd='convert '+filename+' -font ''Arial'' -pointsize 24 -gravity northwest -annotate 0 '+"'"+text1+"'"+' '+save_name
	    #cmd='convert '+filename+' -crop 1582x826+348+488 '+save_name
            cmd='convert '+filename+' -crop 1432x826+408+488 '+save_name

	    os.system(cmd)



end

