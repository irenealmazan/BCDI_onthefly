function [ first ] = extract_number_from_string(string,number,prefix)
%jclark
%extract the number from a string
%retunr the string though (preserves leading zeros)
%but remove any leading whitespace

idx = regexp(char(string),'\d+');
nums = regexp(char(string),'\d+','match');   %get the numbers from the file name

if exist('prefix') == 0
    if exist('number') ~= 0
        first=strtrim(char(nums(number)));        %get first 
    else
        first=strtrim(char(nums(end)));        %get first 
    end
else
    
    first=char(regexp(string,[prefix,'\d+'],'match'));  %nums=char(regexp(char(amp_fs(1).name),'T\d+','match'));
    first(1:numel(prefix))=[];                                  %nums=nums(2:end);
    
end

first=str2num(first); %careful with this
end

