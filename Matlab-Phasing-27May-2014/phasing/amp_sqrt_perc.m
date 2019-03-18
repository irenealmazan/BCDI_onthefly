function [ pns] = amp_sqrt_perc(pn,pow,perc)
%jclark   Detailed explanation goes here

try
    pow;
catch
    pow=0.5;
end
    
try
    perc;
catch
    perc=.02;
end

amp=abs(pn);
ph=atan2(imag(pn),real(pn));
tot=sum(amp(:));

sz=size(amp);

[B IX]=sort(amp(:),'descend');

np=round(perc*prod(sz));

val=amp(IX(np));

amp(amp >= val)=(amp(amp >= val)).^(pow);

%amp=(amp).^(pow);

amp=amp/sum(amp(:))*tot;

pns=amp.*exp(1i*ph);




end

