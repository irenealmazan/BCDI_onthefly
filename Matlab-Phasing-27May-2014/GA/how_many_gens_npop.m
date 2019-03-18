function [ total_npopngen ] = how_many_gens_npop(npop,ngens,nremove,nmin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

npop_gen=zeros([ngens,1]);

for qq=1:ngens
    
    npop_qq=npop-(qq-1)*nremove;
    
    if npop_qq < nmin,npop_qq=npop_gen(qq-1);end
    
    npop_gen(qq)=npop_qq;
    
    disp(num2str(npop_qq))
    
end

total_npopngen=sum(npop_gen(:));

end

