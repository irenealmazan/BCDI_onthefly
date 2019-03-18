function [ pns] = amp_sqrt_th(pn,pow)
%jclark   Detailed explanation goes here

try
    pow;
catch
    pow=0.5;
end

th=.2;

amp=abs(pn);
ph=atan2(imag(pn),real(pn));
tot=sum(amp(:));

%amp=sqrt(amp);
amp=(amp).^(pow);

amp=amp/sum(amp(:))*tot;

amp( abs(pn) < th*max(abs(pn(:))) )=abs(pn( abs(pn) < th*max(abs(pn(:))) ));

pns=amp.*exp(1i*ph);




end

