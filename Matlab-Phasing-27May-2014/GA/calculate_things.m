function [ out ] = calculate_things(pn)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

out.L1_norm=sum(abs(pn(:)));

out.sharpness=sum(abs(pn(:)).^4);

out.summed_phase=sum_phase(pn);

%ph_pn=zero_phase(pn);

SS=shrink_wrap(pn,.2,.5);
out.area=sum(SS(:));
SS=0;




end

