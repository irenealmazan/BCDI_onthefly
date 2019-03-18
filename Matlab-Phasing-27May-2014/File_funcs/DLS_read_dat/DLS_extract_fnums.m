function [first last] = DLS_extract_fnums(fname,prefix)
%extracts .tif file numbers from .dat scan file from Diamond (I-07)

try
    prefix;
catch
    prefix='p100kImage';         %set a default prefix
end

num_format='000000';              %used to determine length of f names      
file_format='.tif';
disp('______________________________________________________________')
disp('Extracting file numbers from .dat DLS scan file....')
disp(['Looking for files of the form - ',prefix,num_format,file_format])

fid=fopen(fname);           %open file
string=fscanf(fid,'%c');    %import as string
fclose(fid);            %close


locs=strfind(string,prefix);

first=string(locs(1):locs(1)+numel(prefix)+numel(file_format)+numel(num_format));
last=string(locs(end):locs(end)+numel(prefix)+numel(file_format)+numel(num_format));
disp('')
disp(['[First,Last] -> [',strtrim(first),',',strtrim(last),']'])
disp('______________________________________________________________')
disp('')

end

