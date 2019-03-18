function [dq_shift_deriv,dq_shift_deriv2] = calc_dq_deriv2_analytical(qbragg,dth_nominal,thBragg)
    % This function calculates the first derivative of the vecotr dq
    % (linking the different positions of the detector at different theta
    % angles (differents points in the rocking curve).

    global  ki_o kf_o

    qbragg_mod = sqrt(qbragg*qbragg');

    dth = dth_nominal ;
    %

    dq_shift_mod_analytic = 2*qbragg_mod*sind(abs(dth/2));%- 2*sind(dth(jj)/2);%



    if dth< 0

        dq_shift_mod_deriv = -qbragg_mod*(pi/180)*cosd(abs(dth/2));

        dq_shift_mod_deriv2 = 0.5*qbragg_mod*(pi/180)^2*sind(abs(dth/2));

        % see summary slide 50:
        dq_shift_x_analytic_unit = -sind(thBragg+dth/2);%sind(thBragg-abs(dth(jj))/2);I don't understand why do I need a minus sign there!
        dq_shift_z_analytic_unit = cosd(thBragg+dth/2);%cosd(thBragg-abs(dth(jj))/2);%2*cosd(thBragg-dth(jj)/2)*sind(-dth(jj)/2)*qbragg_mod;%


        dq_shift_x_analytic_unit_deriv = -(pi/180)*0.5*cosd(thBragg+dth/2);
        dq_shift_z_analytic_unit_deriv = -(pi/180)*0.5*sind(thBragg+dth/2);
        
        dq_shift_x_analytic_unit_deriv2 = (pi/180)^2*0.25*sind(thBragg+dth/2);
        dq_shift_z_analytic_unit_deriv2 = -(pi/180)^2*0.25*cosd(thBragg+dth/2);
        
    else


        dq_shift_mod_deriv = qbragg_mod*(pi/180)*cosd(abs(dth/2));

        dq_shift_mod_deriv2 = -0.5*qbragg_mod*(pi/180)^2*sind(abs(dth/2));

        % see summary slide 50:

        dq_shift_x_analytic_unit = -sind(abs(thBragg)-abs(dth)/2);
        dq_shift_z_analytic_unit = -cosd(abs(thBragg)-abs(dth)/2);%2*cosd(thBragg-dth(jj)/2)*sind(-dth(jj)/2)*qbragg_mod;%
        
        dq_shift_x_analytic_unit_deriv = +(pi/180)*0.5*cosd(abs(thBragg)-abs(dth)/2);
        dq_shift_z_analytic_unit_deriv = -(pi/180)*0.5*sind(abs(thBragg)-abs(dth)/2);
        
        dq_shift_x_analytic_unit_deriv2 = +(pi/180)^2*0.25*sind(abs(thBragg)-abs(dth)/2);
        dq_shift_z_analytic_unit_deriv2 = +(pi/180)^2*0.25*cosd(abs(thBragg)-abs(dth)/2);

    end


    dq_shift_x_deriv = dq_shift_mod_deriv*dq_shift_x_analytic_unit+dq_shift_mod_analytic*dq_shift_x_analytic_unit_deriv;%-2*qbragg_mod*sind(dth(jj)/2)*cosd(thBragg-dth(jj)/2);%sind(thBragg-dth(jj))*qbragg_mod;

    
    dq_shift_z_deriv =  dq_shift_mod_deriv*dq_shift_z_analytic_unit+dq_shift_mod_analytic*dq_shift_z_analytic_unit_deriv;%2*qbragg_mod*sind(dth(jj)/2)*sind(thBragg-dth(jj)/2);%cosd(thBragg-dth(jj))*qbragg_mod;

    dq_shift_deriv = [dq_shift_x_deriv 0 dq_shift_z_deriv];
    
    dq_shift_x_deriv2 = dq_shift_mod_deriv2*dq_shift_x_analytic_unit+dq_shift_mod_deriv*dq_shift_x_analytic_unit_deriv...
        + dq_shift_mod_deriv*dq_shift_x_analytic_unit_deriv + dq_shift_mod_analytic*dq_shift_x_analytic_unit_deriv2;
    
    dq_shift_z_deriv2 = dq_shift_mod_deriv2*dq_shift_z_analytic_unit+dq_shift_mod_deriv*dq_shift_z_analytic_unit_deriv...
        + dq_shift_mod_deriv*dq_shift_z_analytic_unit_deriv + dq_shift_mod_analytic*dq_shift_z_analytic_unit_deriv2;
    
   dq_shift_deriv2 = [dq_shift_x_deriv2 0 dq_shift_z_deriv2];

end