function [labx labz det2] = extract_spec_scan_det2( spec_file,scan_number)
%jclark
textlines=34; %header size

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

det2=temp.data(:,22);

end

function output_string(save_name,string)
%jclark

fid=fopen([save_name],'w');

%fprintf(fid,'\n');
fprintf(fid,'%c',string);
%fprintf(fid,'\n');

fclose(fid);

end

% function newData1=importfile2(fileToRead1)
% %IMPORTFILE2(FILETOREAD1)
% %  Imports data from the specified file
% %  FILETOREAD1:  file to read
% 
% %  Auto-generated by MATLAB on 25-Mar-2013 22:24:32
% 
% DELIMITER = ' ';
% HEADERLINES = 25;
% 
% % Import the file
% newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);
% 
% % Create new variables in the base workspace from those fields.
% vars = fieldnames(newData1);
% for i = 1:length(vars)
%     assignin('base', vars{i}, newData1.(vars{i}));
% end
%     
% end