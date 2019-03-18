function str = load_python_script(script)
% loads a python file as a string for outputting in another location

file=[return_python_script_path(),script];

fid=fopen(file);            %open the file
str=fscanf(fid,'%c');       %read in the file contents

fclose(fid);  



end

