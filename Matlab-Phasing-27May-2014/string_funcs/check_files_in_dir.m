function [nfiles missingf] =check_files_in_dir(ddir)
%jclark

fnames=rdir([ddir,'*.tif']);

nfiles=size(fnames,1);

fnumbs=[];

for qq=1:nfiles
    
    tempname=char(fnames(qq).name);
    
    [dir_out file_out] = extract_dir_from_string( tempname);
    
    fnumbs(qq)=extract_number_from_string(file_out);
    
end



end

