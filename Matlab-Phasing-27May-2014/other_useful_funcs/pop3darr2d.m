function [ arr ] = pop3darr2d( arrin,nz )
%jclark

arr=zeros([size(arrin),nz]);

for qq=1:nz
    arr(:,:,qq)=arrin;
end


end

