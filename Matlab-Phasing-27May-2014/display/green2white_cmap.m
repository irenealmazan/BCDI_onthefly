function [ cmap ] = green2white_cmap(nn)
%jclark

if exist('nn') ~= 1,nn=256;end

cmap=zeros([nn,3]);

cmap(:,2)=1;  %green channel

cmap(:,1)=(1:nn)/nn;

cmap(:,3)=(1:nn)/nn;


end

