% this scripts compares a calculated and an experimental DP:

global ki_o kf_o X Y Z

qbragg = kf_o - ki_o;

% calculated
jj = 27;

Ry = [cosd(-data_exp(jj).dth) 0 sind(-data_exp(jj).dth);
    0 1 0;
    -sind(-data_exp(jj).dth) 0 cosd(-data_exp(jj).dth)];

ki = (Ry * ki_o.').';
kf = (Ry * kf_o.').';

dq_shift_grid = kf - ki - qbragg;

Qterm = exp(1i* dq_shift_grid(1) * X) .* ...
exp(1i* dq_shift_grid(2) * Y) .* ...
exp(1i* dq_shift_grid(3) * Z);

 Psij_1 = sum( rho.*probe.*Qterm, 3);
 Psij = fftshift(fftn(fftshift(Psij_1)));
    
Psij_mod = Psij.*conj(Psij);

figure;
subplot(1,2,1);
imagesc(sqrt(data_exp(jj).I));
axis image;
title(['experimental dp at ' num2str(data_exp(jj).dth)]);

subplot(1,2,2);
imagesc(sqrt(Psij_mod));
axis image;
title(['calculated dp at '  num2str(data_exp(jj).dth_new)]);
        