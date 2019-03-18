function [ output_args ] = batch_change_variables_func(top_dir,script,str_old,str_new)
%jclark
%batch change things inside a file
%top_dir is the dir where to look for the script
%script is the file you want to change
%str_old is a cell array with the old string to search for
%str_new is the cell array with the replacement strings
%eg. 
%batch_change_variables_func([pwd,'/'],'Matlabphasing_ver1_1',{'params.sigma='},{'params.sigma=0.75;%'})

for qq=1:max(size(str_old))
    oldstr=char(str_old(1,qq));
    newstr=char(str_new(1,qq));
    batch_change_variable(top_dir,script,oldstr,newstr);
end



end

