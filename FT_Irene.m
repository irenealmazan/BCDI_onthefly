function [Psi2,Psi2_mod,rock_Psi2] = FT_Irene(data,data_exp,rho)
    % This function performs a FT for each angle of the rocking curve
    
    global ki_o kf_o X Y Z
    
    Psi2 = zeros(size(data));
    Psi2_mod = zeros(size(data));
    rock_Psi2 = zeros(numel(data_exp),1);

    for jj=1:numel(data_exp)
        Ry = [cosd(-data_exp(jj).dth_new) 0 sind(-data_exp(jj).dth_new);
            0 1 0;
            -sind(-data_exp(jj).dth_new) 0 cosd(-data_exp(jj).dth_new)];
        
        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';
        
        qbragg = kf_o-ki_o;
        
        dq_shift = kf - ki - qbragg;
        
        Qterm = exp(1i* dq_shift(1) * X) .* ...
            exp(1i* dq_shift(2) * Y) .* ...
            exp(1i* dq_shift(3) * Z);
        
        
        Psi2(:,:,jj) = fftshift(fftn(fftshift(sum(Qterm.*rho,3))));
        
        Psi2_mod(:,:,jj) = Psi2(:,:,jj).*conj(Psi2(:,:,jj));
        
        rock_Psi2(jj) = sum(sum(sqrt(Psi2_mod(:,:,jj))));
    end


end