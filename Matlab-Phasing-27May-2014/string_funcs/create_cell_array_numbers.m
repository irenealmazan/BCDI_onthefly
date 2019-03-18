function [ nmbs ] = create_cell_array_numbers(nn,strt,nlen)
%jclark

if exist('strt') ~= 1,strt=1;end
if exist('nlen') ~= 1,nlen=[];end

for qq=strt:strt+nn
   
    if isempty(nlen)
        nmbs{qq}=num2str(qq);
    else
        nmbs{qq}=check_numlength(qq,nlen);
    end
end



end

