function [dir_out file_out] = extract_dir_from_string( string)
%jclark


ind=regexp(string,'/');

dir_out=string(1:ind(end));

file_out=string(ind(end)+1:end);

end

