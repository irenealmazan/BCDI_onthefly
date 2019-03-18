function [ str_n ] = load_script_as_text(file,marker)
%jclark

file=strtrim(file);
disp('')
%display('Opening file....')
%disp(file)

try
    fid=fopen(file);            %open the file

    str=fscanf(fid,'%c');       %read in the file contents

    fclose(fid);                %close the file

    %if exist('marker') || isempty(marker) == 1
    if exist('marker') ~= 0
        if isempty(marker) ~= 0
            ind = regexp(str,marker);

            str_n=str(1:ind(1));
        end
    else
        str_n=str;
    end
catch
   disp('Could not open file....') 
   disp(file) 
   str_n=[]; 
end
disp('')   

end

