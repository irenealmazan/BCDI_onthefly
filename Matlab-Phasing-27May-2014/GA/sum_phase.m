function [ summed_ph ] = sum_phase(pn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pn=zero_phase(pn);

ph=atan2(imag(pn),real(pn));

SS=shrink_wrap(abs(pn),.2,.5);

summed_ph=sum( abs(ph(SS > 0)));

end

