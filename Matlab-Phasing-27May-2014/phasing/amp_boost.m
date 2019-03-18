function [ pnb ] = amp_boost(pn,stuff)
%jclark

try
    stuff;
catch
    stuff=stuff_defaults();
end

amp=abs(pn);
ph=atan2(imag(pn),real(pn));

uval=max(amp(:))*stuff.u_val;
lval=max(amp(:))*stuff.l_th;

ind=( amp < uval & amp > lval);

amp(ind)=uval;

pnb=(amp).*exp(i*ph);




end

function stuff=stuff_defaults()

stuff.l_th=0.6; %.6

stuff.u_val=0.8; %.8

end