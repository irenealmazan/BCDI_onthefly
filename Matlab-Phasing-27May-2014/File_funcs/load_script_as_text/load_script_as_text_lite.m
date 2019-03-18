function [ str_n ] = load_script_as_text_lite(file)
%jclark

file=strtrim(file);
disp('')


try
    fid=fopen(file);            %open the file

    str_n=fscanf(fid,'%c');       %read in the file contents

    fclose(fid);                %close the file

catch
   disp('Could not open file....') 
   disp(file) 
   str_n=[]; 
end
disp('')   

end

