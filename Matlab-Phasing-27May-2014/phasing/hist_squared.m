function [ pnh ] = hist_squared(pn,support)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

amp=abs(pn);
totamp=sum(amp(:));
ph=atan2(imag(pn),real(pn));
[h x]=hist(abs(pn(logical(support))));
toth=sum(h(:));
h=h.^2;
h=h/sum(h(:))*toth;

ind=histeq(abs(pn(logical(support))),h);



amp(logical(support))=ind;

amp=amp/sum(amp(:))*totamp;

pnh=amp.*exp(i*ph);


end

