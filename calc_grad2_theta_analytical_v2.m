function grad2_final_theta = calc_grad2_theta_analytical_v2(probe, rho, data, dth_nominal, dth_delta,dthBragg)
    %%% This function calculates the hessian of the error analytically
    
    global  X Y Z ki_o kf_o
    
   qbragg = kf_o - ki_o;
   
   %%% calculate the derivative of dq_shift with respect to theta:
  [dq_shift_deriv,dq_shift_deriv2] = calc_dq_deriv2_analytical(qbragg,dth_nominal,dthBragg);
  
  %%% calcualte the dq_shift
  dth = dth_nominal + dth_delta;
        
  Ry = [cosd(-dth) 0 sind(-dth);
      0 1 0;
      -sind(-dth) 0 cosd(-dth)];
  
  ki = (Ry * ki_o.').';
  kf = (Ry * kf_o.').';
  
  dq_shift = kf - ki - qbragg;
   
  h3 = gcf;
  %display_calc_dqshift_deriv(qbragg,dth_nominal,dq_shift_deriv,dq_shift,dthBragg,h3);
  %display_calc_dqshift_deriv2(qbragg,dth_nominal,dq_shift_deriv2,dq_shift_deriv,dthBragg);
  
  
  Qterm = exp(1i* dq_shift(1) * X) .* ...
      exp(1i* dq_shift(2) * Y) .* ...
      exp(1i* dq_shift(3) * Z);
  
  % intensity of the difracted wave:
  
    Psij = probe.*rho.*Qterm;
    Psij = sum(Psij,3);                     %Radon oper
    Psij = fftshift(fftn(fftshift(Psij)));
    Psij_conj = conj(Psij);
    Psij_mod = Psij.*conj(Psij);


    %%% derivative of Qterm with respect to dq_shift
    deriv_dq_shift_vector = (dq_shift_deriv(1).*X + dq_shift_deriv(3).*Z);
    
    deriv_Qterm_theta = deriv_dq_shift_vector.*Qterm;
    
    deriv2_dq_shift_vector = (dq_shift_deriv2(1).*X + dq_shift_deriv2(3).*Z);
    
    deriv2_Qterm_theta = (deriv2_dq_shift_vector+1i*deriv_dq_shift_vector.^2).*Qterm;
  
    %bmtemp = circshift(probe, round([data(ii).ypos data(ii).xpos data(ii).zpos]/d2_bragg));

   
    % derivative of the diffraction pattern with respect to dq_shift
    Psij_deriv_theta = probe.*rho.*deriv_Qterm_theta;
    
    Psij_deriv_theta = sum(Psij_deriv_theta,3);
    
    Psij_deriv_theta = fftshift(fftn(fftshift(Psij_deriv_theta)));
    
    Psij_deriv_theta_conj = conj(Psij_deriv_theta);
    
    Psij_deriv_theta_mod = Psij_deriv_theta_conj.*Psij_deriv_theta;
    
    Im_Psij_deriv_theta =  imag(Psij_conj.*Psij_deriv_theta);
    
     % second derivative of the diffraction pattern with respect to dq_shift
    Psij_deriv2_theta = probe.*rho.*deriv2_Qterm_theta;
    
    Psij_deriv2_theta = sum(Psij_deriv2_theta,3);
    
    Psij_deriv2_theta = fftshift(fftn(fftshift(Psij_deriv2_theta)));
    
    Im_Psij_deriv2_theta =  imag(Psij_conj.*Psij_deriv2_theta);
    
    %%%% Calculate the gradient (see A. Tripathi et al., Optic express, 22, 1452 (2014):
   
    
    grad_theta_term1 = (sqrt(flipud(data.I))./(Psij.*Psij_conj).^(3/2)).*Im_Psij_deriv_theta.^2;%Im_Psij_deriv_theta./sqrt(Psij_mod);
    
    grad_theta_term2 = (1 -sqrt(flipud(data.I)./(Psij.*Psij_conj))).*Im_Psij_deriv2_theta;
    
    grad_theta = grad_theta_term1 - grad_theta_term2; %%%% NEW 
            
    grad2_final_theta =2*sum(sum(grad_theta))/(numel(grad_theta));
    
    

    %%% check that you calculate the gradient:
    %error_theta_integ = calc_error_theta_singlepos(rho,probe,data,Qterm);
    %display_calc_grad_theta(probe, rho, data, dth_nominal,grad_final_theta,error_theta_integ)
    %display_calc_grad2_theta(probe, rho, data, dth_nominal,grad2_final_theta,grad_final_theta)
     

end