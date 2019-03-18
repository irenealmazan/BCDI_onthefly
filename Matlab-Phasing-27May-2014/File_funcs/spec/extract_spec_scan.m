function [labx labz] = extract_spec_scan( spec_file,scan_number)
%jclark
textlines=34;%33;

%load the spec file
str=load_script_as_text_lite(spec_file);

[dir0 f0]=extract_dir_from_string(spec_file);

%get the scan
ind1=regexp(str,['#S ',num2str(scan_number),' ']);

ind2=regexp(str,['#S ',num2str(scan_number+1),' ']);

str0=str(ind1:ind2);

save_name=[dir0,f0,'-S',num2str(scan_number)];
output_string(save_name,str0);

temp=importfile1(save_name,textlines);

labx=temp.data(:,1);
labz=temp.data(:,2);

end

function output_string(save_name,string)
%jclark

fid=fopen([save_name],'w');

%fprintf(fid,'\n');
fprintf(fid,'%c',string);
%fprintf(fid,'\n');

fclose(fid);

end