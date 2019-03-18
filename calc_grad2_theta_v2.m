function grad2_final_theta = calc_grad2_theta_v2(probe, rho, data, dth_nominal, dth_delta)
    %%% This function calculates the hessian of the error in a very dummy
    %%% way....
    
    
    
    
    
    th_fine_grid = [-0.00005 0.00005];
    
    %%% calculate manually the gradient of dq:
    for jj = 1:numel(th_fine_grid)
        dth_grid(jj) = dth_nominal + dth_delta +  th_fine_grid(jj);
        
       [grad_final_theta(jj)] = calc_grad_theta(probe, rho, data, dth_nominal,th_fine_grid(jj));
        
        
    end
    
    grad2_final_theta = (grad_final_theta(2)-grad_final_theta(1))/(th_fine_grid(2)-th_fine_grid(1));
    
 

end