function ph_support = phase_constraint(pn,phase_range,type)
%generates a support based on a phase range
% for use with phase constrained algorithms
%does Nd

try
    type;
catch
    type='PCr';
end

phase=atan2(imag(pn),real(pn)); %calc the phase

ind1=(phase <= phase_range(2)); %generate a support <= max phase
ind2=(phase >= phase_range(1)); %generate a support >= min phase
ph_support=1e0*ind1.*ind2;       %multiply to obtain both
   

end