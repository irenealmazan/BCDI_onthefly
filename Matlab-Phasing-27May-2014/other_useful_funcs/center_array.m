function [ array xyz] = center_array(array,xyonly)
%jclark
%careful, returns yxz actually

if exist('xyonly') ~= 1,xyonly=0;end

sz=size(array);

nx=sz(1);

dims=ndims(array);

if dims >= 2, ny =sz(2);end
if dims == 3, nz =sz(3);end

dims=ndims(array);

if dims == 3
    %jesse old version
    %[i,j,k]=ind2sub(size(array),find( abs(array) == max(max(max(abs(array))))));
    
    %au edit to use a smooth3 to center
    arnew = smooth3(array,'box',[3 3 3]);
    [i,j,k]=ind2sub(size(arnew),find( abs(arnew) == max(max(max(abs(arnew))))));
    
elseif dims == 2 
    [i,j]=ind2sub(size(array),find( abs(array) == (max(max(abs(array))))));
elseif dims == 1
    [i j]=ind2sub(size(array),find( abs(array) == (max(max(abs(array))))));
end

if mod(nx,2) == 1, xc = (nx+1)/2; else xc =nx/2;end

if dims == 2
    if mod(ny,2) == 1, yc = (ny+1)/2; else yc =ny/2;end
end
if dims == 3
    if mod(ny,2) == 1, yc = (ny+1)/2; else yc =ny/2;end
    if mod(nz,2) ==1, zc = (nz+1)/2; else zc = nz/2;end
end

if dims == 1, array=circshift(array,[-i(1)+xc]);xyz=[-i(1)+xc];end
if dims == 2, array=circshift(array,[-i(1)+xc,-j(1)+yc]);xyz=[-i(1)+xc,-j(1)+yc];end

if xyonly ~= 1,
    if dims == 3, array=circshift(array,[-i(1)+xc,-j(1)+yc,-k(1)+zc]);xyz=[-i(1)+xc,-j(1)+yc,-k(1)+zc];end
else
    if dims == 3, array=circshift(array,[-i(1)+xc,-j(1)+yc,0]);xyz=[-i(1)+xc,-j(1)+yc,-k(1)+zc];end
end


end
%% old version
% function [ array xyz] = center_array( array )
% %UNTITLED11 Summary of this function goes here
% %   Detailed explanation goes here
% 
% sz=size(array);
% 
% nx=sz(1);
% ny=sz(2);
% 
% dims=ndims(array);
% 
% if dims == 3, nz =sz(3);end
% 
% dims=ndims(array);
% 
% if dims == 3
%     [i,j,k]=ind2sub(size(array),find( abs(array) == max(max(max(abs(array))))));
% else 
%     [i,j]=ind2sub(size(array),find( abs(array) == (max(max(abs(array))))));
% end
% 
% if mod(nx,2) == 1, xc = (nx+1)/2; else xc =nx/2;end
% if mod(ny,2) == 1, yc = (ny+1)/2; else yc =ny/2;end
% 
% if dims == 3
%     if mod(nz,2) ==1, zc = (nz+1)/2; else zc = nz/2;end
% end
% 
% if dims == 2, array=circshift(array,[-i(1)+xc,-j(1)+yc]);xyz=[-i(1)+xc,-j(1)+yc];end
% if dims == 3, array=circshift(array,[-i(1)+xc,-j(1)+yc,-k(1)+zc]);xyz=[-i(1)+xc,-j(1)+yc,-k(1)+zc];end
%    
% 
% 
% end
