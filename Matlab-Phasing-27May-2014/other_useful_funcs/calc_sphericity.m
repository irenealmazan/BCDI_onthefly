function [psi psi_n psi_s] = calc_sphericity(array)
%jclark
%calculates sphericity of an array
%by comparing the surface area to volume ratio
%assumes the array is binary

%[ shell ] = calc_shell_3D(array);
shell= bwperim(array);

Ap=sum(shell(:));
Vp=sum(array(:));

psi=(pi).^(1/3)*(6*Vp).^(2/3)/(Ap);

%return a normlaixed version, taking inot account the >1 because of the
%discrete array size
[xx yy zz]=meshgrid(1:size(array,2),1:size(array,1),1:size(array,3));
xx=xx-max(xx(:))/2;yy=yy-max(yy(:))/2;zz=zz-max(zz(:))/2;
rr=sqrt( xx.^2+yy.^2+zz.^2 ); 
rc=(3*Vp/4/pi).^(1/3); %radius if sphere of same volume

disc_sp=(rr <= rc); %the sphere

Vs=sum(disc_sp(:)); %volume sphere
shell_s=bwperim(disc_sp); %get the perimeter
As=sum(shell_s(:)); %get the area

psi_s=(pi).^(1/3)*(6*Vs).^(2/3)/(As);

psi_n=psi/psi_s;

end

