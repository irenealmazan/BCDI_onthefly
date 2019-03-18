function [dq_shift, Psij_mod] = calc_psij_at_dth(dth_new,bmtemp,rho)

    global ki_o kf_o X Y Z
    
    qbragg = kf_o - ki_o;

    Ry = [cosd(-dth_new) 0 sind(-dth_new);
        0 1 0;
        -sind(-dth_new) 0 cosd(-dth_new)];

    ki = (Ry * ki_o.').';
    kf = (Ry * kf_o.').';

    dq_shift = kf - ki - qbragg;

    Qterm = exp(1i* dq_shift(1) * X) .* ...
    exp(1i* dq_shift(2) * Y) .* ...
    exp(1i* dq_shift(3) * Z);

    % calculate the diffraction pattern with the new angles
    Psij = bmtemp.*rho.*Qterm;
    Psij = sum(Psij,3);                     %Radon oper
    Psij = fftshift(fftn(fftshift(Psij)));
    Psij_conj = conj(Psij);

    Psij_mod = Psij.*Psij_conj;
            

end