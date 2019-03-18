function [Pmrho] = InvFT_Irene(data,data_exp,Psi_mod)
    % This function makes the inverse Fourier transform (from the
    % reciprocal space to the real space, for each slice of the imput
    
    global ki_o kf_o X Y Z
    
    Pmrho = zeros(size(data));

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

          Pmrho_dummy = fftshift(ifftn(fftshift(Psi_mod(:,:,jj))));
          Pmrho_dummy2 = repmat(Pmrho_dummy,[1 1 numel(data_exp)])/numel(data_exp);
          Pmrho =  Pmrho + conj(Qterm).*Pmrho_dummy2;


    end
    


end
