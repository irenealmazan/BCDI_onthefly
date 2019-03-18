function [ pna ] =amplitude_constraint(pn,amp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%amp=amp/sum(abs(amp(:)))*sum(abs(pn(:)));

pna=amp.*exp(i*atan2(imag(pn),real(pn)));


end

