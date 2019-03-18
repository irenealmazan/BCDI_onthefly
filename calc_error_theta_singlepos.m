function [err] = calc_error_theta_singlepos(rho,bmtemp,data_exp,Qterm)

    temp = sum( (rho).*bmtemp.*Qterm, 3);
    Psij = fftshift(fftn(fftshift(temp)));  
    
    %%%% NEW ADDITION:
     Psij_mod = Psij.*conj(Psij);
      
    
    err = sqrt(flipud(data_exp.I)) - sqrt(Psij_mod); 
    err = sum(err(:).^2)/numel(err); 
   

end