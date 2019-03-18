function [ pnpc ] = phase_projector(pn,params)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

phi_max=max(params.phase_range);
phi_min=min(params.phase_range);

phase=angle(pn);
amp=abs(pn);

ind=( phase < phi_min);
amp(ind)=cos( (phase(ind)-phi_min) ).*amp(ind);
phase(ind)=phi_min;


ind=( phase > phi_max);

amp(ind)=cos( (phase(ind)-phi_max) ).*amp(ind);
phase(ind)=phi_max;

pnpc=amp.*exp(i*phase);


end

