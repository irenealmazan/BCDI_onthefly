function  output_python_script(save_dir,script)
%will output a copy of the python script, specified by script
%located in the matlab_phasing/python_scripts_for_matlab/ directory

str=load_python_script(script);

file=[strtrim(save_dir),strtrim(script)];

fid=fopen(file,'w');
fprintf(fid,'%c',str);
fclose(fid);

end

