function [strnum] = check_numlength(number,nlen)
%jclark

if ischar(number) ~= 1,strnum=num2str(number);else strnum = number;end

while numel(strnum) < nlen
    
    strnum=['0',strnum];
    
end

end

