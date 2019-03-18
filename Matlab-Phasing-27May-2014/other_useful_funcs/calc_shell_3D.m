function [ shell ] = calc_shell_3D( array ,edgemethod)
%jclark
%used for determining the surface area of a 3d object
%calculates the shell of an array
edgemethod='canny';

sz=size(array,3);

shell=zeros(size(array));

for qq=1:sz
   
   shell(:,:,qq)=edge(array(:,:,qq),edgemethod);
    
end

end

