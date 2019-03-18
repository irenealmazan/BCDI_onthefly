function [ new_array ] = zero_pad_ver3( input,newx,newy,newz )
%jclark

nd=ndims(input);

nn=size(input);

x=newx-nn(2);
y=newy-nn(1);

if ndims(input) == 3,z=newz-nn(3);else z=0;end

nnc=[floor(x/2),ceil(x/2),floor(y/2),ceil(y/2),floor(z/2),ceil(z/2)];

new_array = init_pad(input,nnc);
new_array = init_crop(new_array,nnc);

end

