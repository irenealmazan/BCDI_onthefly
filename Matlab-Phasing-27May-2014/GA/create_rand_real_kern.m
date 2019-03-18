function [kern Fkern]= create_rand_real_kern(kn,nn,custom)

if exist('custom')
    A=reshape(custom,[kn,kn]); 
else
    A=random('uniform',-1,1,[kn,kn]);
end

A=zero_pad_ver3(A,nn,nn);

F=fftshift(fftn(fftshift(A)));

kern=fftshift(ifftn(fftshift(abs(F).^2)));

Fkern=abs(F).^2;

end