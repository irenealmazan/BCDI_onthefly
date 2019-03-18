%%% This function displays the calculated gradient and check that this is the
%%% slope of the err function vs theta

function [dth_grid,err_thin_grid,tangent_grad] = display_calc_grad_theta(probe, rho, data, dth_nominal,grad_calc,err_0)


    global  X Y Z ki_o kf_o
    
    qbragg = kf_o - ki_o;
    
    th_fine_grid = [-0.05:0.5e-3 :0.05];%[-0.005:1e-5 :0.005];

    for jj = 1:numel(th_fine_grid)
        dth_grid(jj) = dth_nominal  +  th_fine_grid(jj);

        Ry = [cosd(-dth_grid(jj)) 0 sind(-dth_grid(jj));
            0 1 0;
            -sind(-dth_grid(jj)) 0 cosd(-dth_grid(jj))];

        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';

        dq_shift_grid(jj,:) = kf - ki - qbragg;

        Qterm = exp(1i* dq_shift_grid(jj,1) * X) .* ...
        exp(1i* dq_shift_grid(jj,2) * Y) .* ...
        exp(1i* dq_shift_grid(jj,3) * Z);
        
        err_thin_grid(jj) = calc_error_theta_singlepos(rho,probe,data,Qterm);

    end
    
       %figure;
       %h1 = plot(dth_grid,err_thin_grid);
       %hold on;
       cstx = err_0-grad_calc*(dth_nominal);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
       tangent_grad = cstx + grad_calc* dth_grid;
       
       
     

end