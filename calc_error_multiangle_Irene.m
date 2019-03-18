function [errtot] = calc_error_multiangle_Irene(probe, rho, data)

global d2_bragg X Y Z

errtot=0;

for ii=1:numel(data)
    
    Qterm = exp(i* data(ii).dqshift(1) * X) .* ...
            exp(i* data(ii).dqshift(2) * Y) .* ...
            exp(i* data(ii).dqshift(3) * Z);
        
    %bmtemp(ii).A = circshift(probe, round([data(ii).ypos data(ii).xpos data(ii).zpos]/d2_bragg));
   

    temp = sum( rho.*probe(ii).A.*Qterm, 3);
    Psij = fftshift(fftn(fftshift(temp)));
    
    Psij_mod = Psij.*conj(Psij);
    
        
    err = sqrt(data(ii).I) - sqrt(Psij_mod);
    err = sum(err(:).^2)/numel(err);        
    errtot = errtot + err;    
end
