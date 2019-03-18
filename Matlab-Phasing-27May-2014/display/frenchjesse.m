function [ cmap ] = frenchjesse()
%jclark

%if exist('nn') ~= 1,nn=256;end

cmap=zeros([256,3]);

%r 1 -> 0 -> 0
%g 1 -> 1 -> 0
%b 1 -> 0 -> 0

v1=85;
v2=150;
v3=236;

ind1=(1:v1);
ind2=((v1+1):v2);
ind3=((v2+1):v3);

cmap(ind1,1)=reverse((1:numel(ind1))/numel(ind1));

cmap(ind1,2)=1;
cmap(ind2,2)=1;
cmap(ind3,2)=reverse((1:numel(ind3))/numel(ind3));

cmap(ind1,3)=reverse((1:numel(ind1))/numel(ind1));



end