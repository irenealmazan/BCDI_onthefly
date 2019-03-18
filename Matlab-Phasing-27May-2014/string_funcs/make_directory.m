function make_directory(additional_dirs)
%jclark
%make additional dirsctorys from a base dir
%
%adidtional_dirs is a cell array



for qq=1:size(additional_dirs,2)
    
    temp=[char(additional_dirs(qq))];
    if isdir(temp) ~= 1,mkdir(temp);end
    
end



end

