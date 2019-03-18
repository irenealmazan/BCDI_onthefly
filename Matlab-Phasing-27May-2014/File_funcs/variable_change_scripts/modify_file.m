function  modify_file(file,oldstr,newstr )
%this will load the text file in file and replace the string
%in oldstr with whats in newstr.  useful for changing variables
%in matlab_phasing without having to open it.  good for doing batch jobs
file=strtrim(file);
disp('')
display('Opening file....')
disp(file)
try
    fid=fopen(file);            %open the file

    str=fscanf(fid,'%c');       %read in the file contents

    fclose(fid);                %close the file

    disp(['Replacing ''',oldstr,''' with ''',newstr,''])  

    mstr = strrep(str,oldstr,newstr); %replace string

    fid=fopen(file,'w');

    fprintf(fid,'%c',mstr);

    fclose(fid);
catch
   disp('Could not open file....') 
   disp(file) 
    
end
disp('')    
end

