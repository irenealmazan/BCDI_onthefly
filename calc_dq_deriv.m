function [dq_shift_deriv] = calc_dq_deriv(dth_nominal,dth_delta)
    % This function calculates the first derivative of the vecotr dq
    % (linking the different positions of the detector at different theta
    % angles (differents points in the rocking curve). 

   global  ki_o kf_o
    
    qbragg = kf_o - ki_o;
    
   th_fine_grid = [-0.00005 0.00005];
   
   %%% calculate manually the gradient of dq:
   for jj = 1:numel(th_fine_grid)
        dth_grid(jj) = dth_nominal + dth_delta +  th_fine_grid(jj);
        
        Ry = [cosd(-dth_grid(jj)) 0 sind(-dth_grid(jj));
            0 1 0;
            -sind(-dth_grid(jj)) 0 cosd(-dth_grid(jj))];
        
        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';
                
        dq_shift_grid(jj,:) = kf - ki - qbragg;
        
       
   end
    
  dq_shift_deriv(1) = (dq_shift_grid(2,1)-dq_shift_grid(1,1))/(th_fine_grid(2)-th_fine_grid(1));
           
  dq_shift_deriv(2) = (dq_shift_grid(2,2)-dq_shift_grid(1,2))/(th_fine_grid(2)-th_fine_grid(1));
  
  dq_shift_deriv(3) = (dq_shift_grid(2,3)-dq_shift_grid(1,3))/(th_fine_grid(2)-th_fine_grid(1));

  %%% check the position of dq_shift_grid
  %{
      Qterm_up = exp(1i* dq_shift_grid(1,1) * X) .* ...
          exp(1i* dq_shift_grid(1,2) * Y) .* ...
          exp(1i* dq_shift_grid(1,3) * Z);

      error_theta_integ_up = calc_error_theta_singlepos(rho,probe,data,Qterm_up);

      Qterm_down = exp(1i* dq_shift_grid(2,1) * X) .* ...
          exp(1i* dq_shift_grid(2,2) * Y) .* ...
          exp(1i* dq_shift_grid(2,3) * Z);

      error_theta_integ_down = calc_error_theta_singlepos(rho,probe,data,Qterm_down);
  
      hold on; plot(dth_grid(1),error_theta_integ_up,'or'); plot(dth_grid(2),error_theta_integ_down,'or');
  %}
  
 %%% check graphically that this the derivative of dq/dtheta:
  %{ 
  th_fine_grid = [-0.1:1e-2 :0.1];
   
   for jj = 1:numel(th_fine_grid)
        dth_grid(jj) = dth_nominal + dth_delta +  th_fine_grid(jj);
        
        Ry = [cosd(-dth_grid(jj)) 0 sind(-dth_grid(jj));
            0 1 0;
            -sind(-dth_grid(jj)) 0 cosd(-dth_grid(jj))];
        
        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';
                
        dq_shift_grid(jj,:) = kf - ki - qbragg;
        
       
   end
   
   figure;subplot(121);plot(dth_grid,dq_shift_grid(:,1));
   hold on;      
   cstx = dq_shift(1)-dq_shift_deriv(1)*(dth_nominal + dth_delta);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
   plot(dth_grid,cstx + dq_shift_deriv(1)* dth_grid);


   subplot(122);plot(dth_grid,dq_shift_grid(:,3));
   hold on;      
   cstx = dq_shift(3)-dq_shift_deriv(3)*(dth_nominal + dth_delta);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
   plot(dth_grid,cstx + dq_shift_deriv(3)* dth_grid);
  %}
  
    %%% test of the derivative with respect to theta: calculate
    %%% analytically the derivative dq/dtheta following Eqs. in the report
    %{
    for jj = 1:numel(dth_delta)
        dth(jj) = dth_nominal + dth_delta(jj);
        
        Ry = [cosd(-dth(jj)) 0 sind(-dth(jj));
            0 1 0;
            -sind(-dth(jj)) 0 cosd(-dth(jj))];
        
        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';
                
        dq_shift(jj,:) = kf - ki - qbragg;
        
        dq_shift_x_analytic(jj) = sind(dth(jj)/2)*sind(thBragg-dth(jj)/2)*qbragg_mod; % taken from Eq. 28 of pos_grad_v1
        
        dq_shift_z_analytic(jj) = -sind(dth(jj)/2)*cosd(thBragg-dth_delta(jj)/2)*qbragg_mod;
                
        dq_shift_x_deriv(jj) = sind(thBragg-dth(jj))*qbragg_mod;
        
        dq_shift_z_deriv(jj) = -cosd(thBragg)*qbragg_mod;
        
        Qterm_x(jj) = exp(1i* dq_shift(jj,1) * X(1,1,1));
        
        deriv_Qterm_x_theta(jj) = 1i*(sind(dth_nominal-dth_delta(jj)).* X(1,1,1)).*qbragg_mod.*Qterm_x(jj);
        
    end
    %}




end