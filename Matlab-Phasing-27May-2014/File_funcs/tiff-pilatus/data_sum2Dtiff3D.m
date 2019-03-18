%use to specify file numbers to sum into 3d tiff's
%either specify the complete file names ie.
% files={first1,last1,first2,last2,.....}
%e.g.
%files={'Au0311_09141.tif','Au0311_09241.tif','Au0311_09243.tif','Au0311_09323.tif'};

%or the file prefix,suffix and numbers.  be careful with leading 0's. (i.e.
% put them in the prefix, see below)

dat_f_nums=[9141,9241,9243,9323];
files=create_fname_numb('Au0311_0','.tif',dat_f_nums);
%%
name_of_this_file='data_sum2Dtiff3D';
dir_file=which(name_of_this_file);    %finds the current folder 
this_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
save_dir=this_dir;
accum=1;



nfiles=size(files);             %get the number of files
nfiles=max(nfiles);
sets=nfiles/2;

%for qq = 1:sets

for qq = 1:2:nfiles
    
    disp('summing....')
    one=char(files(qq));
    two=char(files(qq+1));
    file_load={one,two};
    Tiff2D_2_Tiff3D( file_load,accum,save_dir)
    
    
end