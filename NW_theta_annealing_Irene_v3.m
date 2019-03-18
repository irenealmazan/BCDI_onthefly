function [dth_new,grad_final_theta,dq_shift,errlist_theta] = NW_theta_annealing_Irene_v3(probe, rho,data_exp,Niter_theta,index_to_distort,dthBragg,cnt_theta)
    %%% In this function we constrain the shift in theta to be the same for
    %%% all the points in the grid. 
    
    global X Y Z ki_o kf_o d2_bragg

    qbragg = kf_o - ki_o;
    
    %{
    % This is the last bit of the script NW_calc_rocking_curve.m

    thscangrid_rock = [-.01:.001:.01];

    thscanlist_rock = zeros(numel(thscangrid_rock)*numel(thscanlist), 1);

    dqlist_rock = zeros(numel(thscangrid_rock)*numel(thscanlist), 3);


    cnt = 1;

    for ii=1:numel(thscanlist)

        for jj=1:numel(thscangrid_rock)

            thscan = thscanvals(ii) +  thscangrid_rock(jj);

            thscanlist_rock(cnt) = thscan;

            Ry = [cosd(thscan) 0 sind(thscan);
                0 1 0;
                -sind(thscan) 0 cosd(thscan)];
            ki = (Ry *ki_o')';
            kf = (Ry *kf_o')';

            dqlist_rock(cnt,:) = (kf-ki) - qbragg;

            cnt = cnt + 1;
        end

    end

    %}

   % figure(401);clf;setfigsize(gcf, 1200,500);
   
    dq_shift = zeros(numel(data_exp),3);
   
    
    dth_new = zeros(numel(data_exp),1);
    for ii=1:numel(data_exp)
         dth_new(ii) = data_exp(ii).dth_new;
    end
    
    
     dir_name = ['results/trial8/iter' num2str(cnt_theta)];
     system(sprintf('mkdir %s',dir_name));
    
    errlist_theta = zeros(numel(data_exp),Niter_theta);
          
    for ntheta = 1:Niter_theta
        
       %orderrandom = randperm(numel(data_exp));
        
       %orderrandom = randperm(numel(index_to_distort));
       
        for ii = index_to_distort%1:numel(data_exp)%index_to_distort(orderrandom)%
            
           
            %bmtemp = circshift(probe, round([data_exp(ii).yshift_delta data_exp(ii).xshift_delta data_exp(ii).zshift_delta]/d2_bragg));
            
            bmtemp = probe(ii).A;
            
            display(['dth_new = ' num2str(dth_new(ii))])
            
            % numerical calculation of the gradient
            [grad_final_theta,err] = calc_grad_theta(bmtemp, rho, data_exp(ii), dth_new(ii),0,dthBragg);
            
            % numerical calculation of the Hessian:
            [grad2_final_theta_num] = calc_grad2_theta_v2(bmtemp, rho, data_exp(ii), dth_new(ii),0);
            
            % calculate the Hessian (see Nocedal Eq. 3.2 page 35)
            %[grad2_final_theta_analytic] = calc_grad2_theta_analytical_v2(bmtemp, rho, data_exp(ii), dth_new(ii),0,dthBragg);
            
            h1 = figure;
            
            % check that the gradients are well calculated
            [dth_grid,err_thin_grid,tangent_grad] = display_calc_grad_theta(bmtemp, rho, data_exp(ii), dth_new(ii),grad_final_theta,err);
            
            subplot(151); 
            plot(dth_grid,err_thin_grid);
            hold on;
            plot(dth_grid,tangent_grad,'r');
            title(['error around ' num2str(data_exp(ii).dth)]);
            
            [dth_grid,grad_final_theta_thingrid,tangent_grad2] = display_calc_grad2_theta(bmtemp, rho, data_exp(ii), dth_new(ii),grad2_final_theta_num,grad_final_theta,dthBragg);
            subplot(152); 
            plot(dth_grid,grad_final_theta_thingrid);
            hold on;
            plot(dth_grid,tangent_grad2,'r');
            title(['first derivative of error around ' num2str(data_exp(ii).dth)]);
           
            
            % see Nocedal Eq. 3.2 page 35 This is the Newton approach for the
            % steepest descent:
            theta_step = -2*grad_final_theta/grad2_final_theta_num;
            
            [dq_shift_before, Psij_mod_before] = calc_psij_at_dth(dth_new(ii),bmtemp,rho);
            
            subplot(153);
            imagesc(sqrt(Psij_mod_before));
            axis image;
            title('calculated dp with non corrected angles');
            
            dth_new(ii) = theta_step + dth_new(ii);
            
            [dq_shift(ii,:), Psij_mod_after] = calc_psij_at_dth(dth_new(ii),bmtemp,rho);
            
            errlist_theta(ii,ntheta) = err;
            % calculat the new data_exp(ii).dqshift
            
           
            subplot(154);
            imagesc(sqrt(Psij_mod_after));
            axis image;
            title('calculated dp with corrected angles');
            
            subplot(155);
            imagesc(sqrt(data_exp(ii).I));
            axis image;
            title('experimental dp ');
            
             
            %save the figures to track the evolution of the error metric
            %with iterations:
           
            savefig(h1,[dir_name '/angle' num2str(ii)]);
           
            
        end
        
    end
    
        
        %{
    %     Qterm = exp(i* dqlist_rock(dthsearchind(ind),1) * X) .* ...
    %             exp(i* dqlist_rock(dthsearchind(ind),2) * Y) .* ...
    %             exp(i* dqlist_rock(dthsearchind(ind),3) * Z);
    %     temp = sum( abs(NW).*bmtemp.*Qterm, 3);
    %     Psij = fftshift(fftn(fftshift(temp)));
    %
    %     subplot(132); imagesc(abs(Psij)); axis image;
    %     subplot(133); imagesc(sqrt(data_exp(ii).I)); axis image
    %     drawnow;
    %
    %     data_exp(ii).dqx_fit = dqlist_rock(dthsearchind(ind),1);
    %     data_exp(ii).dqy_fit = dqlist_rock(dthsearchind(ind),2);
    %     data_exp(ii).dqz_fit = dqlist_rock(dthsearchind(ind),3);
    %     data_exp(ii).dth_fit = thscanlist_rock(dthsearchind(ind));
    %
    %     thfit_chan = [thfit_chan thscanlist_rock(dthsearchind(ind))];
        %}
    
      %%% TEST OF  GRADIENT
        %{
            for jj= dthsearchind(2:end-1)
                grad_manual = test_grad_theta_manually(jj,thscan,fitQerr,data_exp(ii),bmtemp,rho,grad_final_theta(ii).grad(jj));

            end
        
        %}
    
    
end