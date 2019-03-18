function output_string(save_name,string)
%jclark

fid=fopen([save_name],'w');

fprintf(fid,'\n');
fprintf(fid,'%c',string);
fprintf(fid,'\n');

fclose(fid)

end
