function [ output_args ] = imagemagick_annotate( input_args )
%jclark
%use imagemagick to crop images

%recursive string
string_f='*.png';

%get current dir
name_of_this_file='imagemagick_annotate';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
fdir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

%search for the file names
fnames=rdir([fdir,'**/*',string_f]);

for qq=1:size(fnames,1)
         
    filename=fnames(qq).name; 
 	
    [fdir fname]=extract_dir_from_string(filename);
    
    save_name=[fdir,'A-',fname];
 	    
    disp(filename)
 	    
    str=['''',num2str(round(qq*20)),' mins'' " '];
    
    cmd=['convert ',filename,' -font ''Arial'' -pointsize 30 -fill black -draw "text 10,30 ',str,save_name];

    %convert 230.png -font 'Arial' -pointsize 120 -fill black -draw "text 50,150 'x-y'" -draw "text 1600,150 'x-z'" -draw "text 2900,150 'z-y'" -draw "text 50 1000 '+230 ps'" 230t.png

    
 	system(cmd)
 
end

end

