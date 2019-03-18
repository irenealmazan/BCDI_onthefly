function batch_change_variable(dir,file,oldstr,newstr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

files=rdir([[dir],'**/*',file]);

fnames=char(files.name);

nn=size(fnames,1);

for qq = 1:nn
    modify_file(fnames(qq,:),oldstr,newstr );
end

end

