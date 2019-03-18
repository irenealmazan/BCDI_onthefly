function [ pns] = amp_sqrt_val(pn,pow,val)
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

val=val*max(amp(:));

%amp(amp >= val)=(amp(amp >= val)).^(pow);
amp(amp >= val)=mean(mean(mean(amp(amp >= val))));
%amp=(amp).^(pow);

amp=amp/sum(amp(:))*tot;

pns=amp.*exp(1i*ph);




end

