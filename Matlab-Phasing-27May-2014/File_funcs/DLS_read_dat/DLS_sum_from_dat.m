function DLS_sum_from_dat
%creates 3d tiff stacks by reading the .dat file from the scan output
%produced at Diamond Light Source.  Just enter as many .dat files as
%you like and it will search for the .tif file names then create a 3D tiff
%stack (32bit) for each .dat file.
%change the prefix to the constant part of the filename before the number
%if the .tif file names are known (first and last) then use
% data_sum2Dtiff3D.m instead.
%%
%1.
dat_f_nums=[30769,30779,30787,30795,30798,30808,30811,30823,30834,30846,30849, ...
    30852,30855,30858,30861,30864,30867,30873];
%2.
dat_files=create_fname_numb('','.dat',dat_f_nums);

%can alternatively just write the files explcitily (see line below)
%dat_files={'30703.dat'};            %.dat file to extract file numbers from

%3.
prefix='p100kImage';    %prefix for .tif files to load
save_prefix='G42-C3-';                     %can specify a prefix for the output files
accum=1;                            %accumulations per theta.  
%% probably won't need to change anything below here

name_of_this_file='DLS_sum_from_dat';
dir_file=which(name_of_this_file);    %finds the current folder 
this_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory    
save_dir=this_dir;

save_dir=[save_dir,save_prefix];

for qq = 1:max(size(dat_files))
    
    [first last] = DLS_extract_fnums(char(dat_files{qq}),prefix);   %get .tif file numbers
    first=strtrim(first);
    last=strtrim(last);
    
    aa=rdir(['**/*',first]);
    bb=rdir(['**/*',last]);
    first=aa.name;
    last=bb.name;
    
    if qq == 1,files={first,last};else files=[files,{first,last}];end
    
end

nfiles=size(files);             %get the number of files
nfiles=max(nfiles);
sets=nfiles/2;
counter=0;

for qq = 1:2:nfiles
    counter=counter+1;
    disp('summing....')
    one=char(files(qq));
    two=char(files(qq+1));
    file_load={one,two};
    Tiff2D_2_Tiff3D( file_load,accum,[save_dir,'#',num2str(dat_f_nums(counter)),'-'])
    
    
end



end

function names=create_fname_numb(prefix,suffix,numbers)

nfiles=max(size(numbers));

names='';

for qq = 1:nfiles
    name=[prefix,num2str(numbers(qq)),suffix];
    names=[names,{name}];
    
end




end


