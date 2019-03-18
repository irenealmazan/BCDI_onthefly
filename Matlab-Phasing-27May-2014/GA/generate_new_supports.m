function [ new_s ] = generate_new_supports(pn_g,threshold,sigma,type)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

try
    type;
catch
    type='gauss';
end

new_s=zeros([size(pn_g)]);

if ndims(pn_g) == 3,npop=size(pn_g,3);end
if ndims(pn_g) == 4,npop=size(pn_g,4);end

for qq=1:npop
    
   if ndims(pn_g) == 3,
       new_s(:,:,qq)=shrink_wrap(pn_g(:,:,qq),threshold,sigma,type);
   end
       
   if ndims(pn_g) == 4,
       new_s(:,:,:,qq)=shrink_wrap(pn_g(:,:,:,qq),threshold,sigma,type);
   end
       
    
end



end

