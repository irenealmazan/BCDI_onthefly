function [grad_final_theta,error_theta_integ] = calc_grad_theta(probe, rho, data, dth_nominal, dth_delta,thBragg)
    % this function calculates the gradient of the error metric with
    % respect to theta. dth_delta = 0 if we are using the function to
    % refine the angular position. 

   global  X Y Z ki_o kf_o
    
   qbragg = kf_o - ki_o;
   
   %%% calculate the derivative of dq_shift with respect to theta:
  %[dq_shift_deriv] = calc_dq_deriv2_analytical(qbragg,dth_nominal,thBragg);
   [dq_shift_deriv] = calc_dq_deriv(dth_nominal,0);
  
  %%% calcualte the dq_shift
  dth = dth_nominal + dth_delta;
        
  Ry = [cosd(-dth) 0 sind(-dth);
      0 1 0;
      -sind(-dth) 0 cosd(-dth)];
  
  ki = (Ry * ki_o.').';
  kf = (Ry * kf_o.').';
  
  dq_shift = kf - ki - qbragg;
   
  %h2 = figure;
  %display_calc_dqshift_deriv(qbragg,dth_nominal,dq_shift_deriv,dq_shift,thBragg,h2)

  
  Qterm = exp(1i* dq_shift(1) * X) .* ...
      exp(1i* dq_shift(2) * Y) .* ...
      exp(1i* dq_shift(3) * Z);
  

    %%% derivative of Qterm with respect to dq_shift
    deriv_Qterm_theta = (dq_shift_deriv(1).*X + dq_shift_deriv(3).*Z).*Qterm;
  
    %bmtemp = circshift(probe, round([data(ii).ypos data(ii).xpos data(ii).zpos]/d2_bragg));

    Psij = probe.*rho.*Qterm;
    Psij = sum(Psij,3);                     %Radon oper
    Psij = fftshift(fftn(fftshift(Psij)));
    Psij_conj = conj(Psij);
    

   
    % derivative of the diffraction pattern with respect to dq_shift
    Psij_deriv_theta = probe.*rho.*deriv_Qterm_theta;
    
    Psij_deriv_theta = sum(Psij_deriv_theta,3);
    
    Psij_deriv_theta = fftshift(fftn(fftshift(Psij_deriv_theta)));
    
    Im_Psij_deriv_theta = imag(Psij_conj.*Psij_deriv_theta);
    
    %%%% Calculate the gradient (see A. Tripathi et al., Optic express, 22, 1452 (2014):
     Psij_mod = Psij.*conj(Psij);
    
    
    grad_theta = (1 -sqrt(flipud(data.I)./(Psij.*Psij_conj))).* Im_Psij_deriv_theta; %%%% NEW 
            
    grad_final_theta = -2*sum(sum(grad_theta))/(numel(grad_theta));
    
    

    %%% check that you calculate the gradient:
    error_theta_integ = calc_error_theta_singlepos(rho,probe,data,Qterm);
    %display_calc_grad_theta(probe, rho, data, dth_nominal,grad_final_theta,error_theta_integ)
    %display_calc_grad2_theta(probe, rho, data, dth_nominal,grad2_final_theta,grad_final_theta)
    
end