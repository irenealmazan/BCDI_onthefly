function dir_py = return_analyse_temp_series_script_path
%returns the path of where this file is located

name_of_this_file='return_analyse_temp_series_script_path';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir_py=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory


end

