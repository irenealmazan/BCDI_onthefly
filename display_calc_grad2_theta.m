function [dth_grid,grad_final_theta,tangent_grad2] = display_calc_grad2_theta(probe, rho, data, dth_nominal ,grad2_calc,grad_calc_0,dthBragg)
%%% This function displays the calculated Hessian second deriv) and check that this is the
%%% slope of the gradient function vs theta


    global  X Y Z ki_o kf_o

    qbragg = kf_o - ki_o;

    th_fine_grid = [-0.05:1e-3 :0.05];

    for jj = 1:numel(th_fine_grid)

        dth_grid(jj) = dth_nominal  +  th_fine_grid(jj);

        [grad_final_theta(jj)] = calc_grad_theta(probe, rho, data, dth_nominal  +  th_fine_grid(jj),0,dthBragg);
    end

    %figure;
    %h2 = plot(dth_grid,grad_final_theta);
    %hold on;
    cstx = grad_calc_0-grad2_calc*(dth_nominal);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
    tangent_grad2 = cstx + grad2_calc* dth_grid;
    %plot(dth_grid,cstx + grad2_calc* dth_grid);




end