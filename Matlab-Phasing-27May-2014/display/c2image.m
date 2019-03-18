function im = c2image(a, varargin)
% FUNCTION IM = C2IMAGE(A)
% 
% Returns a RGB image of complex array A where
% the phase is mapped to hue, and the amplitude
% is mapped to brightness.

absa = abs(a);
phasea = angle(a);

% (optional second argument can switch between various plotting modes)
abs_range = [];
if nargin==2
    m = varargin{1};
elseif nargin==3
    m = varargin{1};
    abs_range = varargin{2};
else
    m = 1;
end

if isempty(abs_range)
    nabsa = absa/max((absa(:)));
else
    nabsa = (absa - abs_range(1))/(abs_range(2) - abs_range(1));
    nabsa(nabsa < 0) = 0;
    nabsa(nabsa > 1) = 1;
end

switch m
  case 1
    im_hsv = zeros([size(a) 3]);
      
    if ndims(a) == 2
        im_hsv(:,:,1) = mod(phasea,2*pi)/(2*pi);
        im_hsv(:,:,2) = 1;
        im_hsv(:,:,3) = nabsa;
        im = hsv2rgb(im_hsv);
    else
        
        im_hsv(:,:,:,1) = mod(phasea,2*pi)/(2*pi);
        im_hsv(:,:,:,2) = 1;
        im_hsv(:,:,:,3) = nabsa;
        im = hsv2rgb(im_hsv);
  
    end
    
  case 2
      if ndims(a) == 2
        im_hsv = ones([size(a) 3]);
        im_hsv(:,:,1) = mod(phasea,2*pi)/(2*pi);
        im_hsv(:,:,2) = nabsa;
        im = hsv2rgb(im_hsv);
      else
        im_hsv = ones([size(a) 3]);
        im_hsv(:,:,:,1) = mod(phasea,2*pi)/(2*pi);
        im_hsv(:,:,:,2) = nabsa;
        im = hsv2rgb(im_hsv);
      end
end

