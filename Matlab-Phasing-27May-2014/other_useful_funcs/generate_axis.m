function [ axis ] = generate_axis(nx,ny,nz,ax)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

axis=zeros([ny,nx,nz]);

sz2=round([ny nx nz]/2);

axis(sz2(1):sz2(1)+ax,sz2(2),sz2(3))=2;
axis(sz2(1),sz2(2):sz2(2)+ax,sz2(3))=1;
axis(sz2(1),sz2(2),sz2(3):sz2(3)+ax)=3;

end

