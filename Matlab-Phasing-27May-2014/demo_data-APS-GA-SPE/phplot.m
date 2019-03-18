function varargout=phplot(varargin)

%PHPLOT(FIELD)
%Plots the phase of FIELD in a continuous color scale (hue) and represents
%the normalized amplitude as brightness (r+g+b)*amplitude.
%PHPLOT(FIELD,AMP,FLAG)
%If AMP = 0 the amplitude is not plot
%If FLAG = 1 the function creates a figure with a dial scale (from 0 to
%2*pi) and radial brightness (from 0 to one)
%A=PHPLOT(...) creates a 3D uint8 array that can be saved as an image with
%IMWRITE(A,'filename','fmt').
%Iacopo Mochi, Lawrence Berkeley National Laboratory 06/6/2010

switch nargin
    case 1
        field=varargin{1};
        Amp=1;
        scale=0;
    case 2
        field=varargin{1};
        Amp=varargin{2};
        scale=0;
    case 3
        field=varargin{1};
        Amp=varargin{2};
        scale=varargin{3};
    case 0
        print('PHPLOT requires at least 1 input argument')
        exit
end

Im=imag(field);
Re=real(field);

phase=atan2(Im,Re);
amplitude=abs(field);
amplitude=amplitude/max(amplitude(:));
if Amp==0
    amplitude=ones(size(amplitude));
end
A=zeros(size(field,1),size(field,2),3);     %Declare RGB array

A(:,:,1)=0.5*(sin(phase)+1).*amplitude;     %Red
A(:,:,2)=0.5*(sin(phase+pi/2)+1).*amplitude;%Green
A(:,:,3)=0.5*(-sin(phase)+1).*amplitude;    %Blue

image(A)

A=uint8(A*255);

if scale>0
    figure
    phase=-pi:2*pi/255:pi;
    r=0.5*(sin(phase)+1);     %Red
    g=0.5*(sin(phase+pi/2)+1);%Green
    b=0.5*(-sin(phase)+1);    %Blue
    
    warphase=[r(:)';g(:)';b(:)'];
    
    [x,y]=meshgrid(-1:2/255:1);
    a=(1i*y+x);
    
    colormap(warphase');
    a((x.^2+y.^2)>1)=nan;
    phplot(a);
    axis image
    axis off
    
end

switch nargout
    case 1
        varargout{1}=A;
    otherwise
end

%example

% [x,y]=meshgrid(-3*pi:0.1:3*pi); 
% F=sin(x.*y)-x*i; 
% subplot(1,3,1) 
% imagesc(abs(F)) 
% subplot(1,3,2) 
% imagesc(angle(F)) 
% subplot(1,3,3) 
% phplot(F,1,0); 
% axis image