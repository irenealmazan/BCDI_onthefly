function [ pns] = amp_sqrt(pn,pow)
%jclark   Detailed explanation goes here

try
    pow;
catch
    pow=0.5;
end
    

amp=abs(pn);
ph=atan2(imag(pn),real(pn));
tot=sum(amp(:));

%amp=sqrt(amp);
amp=(amp).^(pow);

amp=amp/sum(amp(:))*tot;

pns=amp.*exp(1i*ph);




end

