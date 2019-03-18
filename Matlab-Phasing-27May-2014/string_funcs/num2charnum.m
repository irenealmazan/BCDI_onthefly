function [ num ] = num2charnum( numarr )
%jclark
%takes an array with numbers (stored as numbers)
%and turns it into a string e.g [1,5,34] -> '1-5-34'

n_ele=max(size(numarr));
num=[];

for qq=1:n_ele
    
   num=[num,num2str(numarr(qq)),'-'];
    
end

num(end)=[];

end

