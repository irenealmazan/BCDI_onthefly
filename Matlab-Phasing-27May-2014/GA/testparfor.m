function [ output_args ] = testparfor( input_args )
%jclark
matlabpool(6)

p=random('uniform',0,1,[64 64 1000]);

tic
for qq=1:1000

    c(:,:,qq)=fftn(p(:,:,qq));
    
end
toc

matlabpool close

end

