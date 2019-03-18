function gradtot = calc_grad_singleangle_Irene(probe, rho, data)

global X Y Z


    Qterm = exp(i* data.dqshift(1) * X) .* ...
        exp(i* data.dqshift(2) * Y) .* ...
        exp(i* data.dqshift(3) * Z);

    bmtemp = probe;
    %bmtemp = circshift(probe, round([data(ii).ypos data(ii).xpos data(ii).zpos]/d2_bragg));

    Psij = bmtemp.*rho.*Qterm;
    Psij = sum(Psij,3);                     %Radon oper
    Psij = fftshift(fftn(fftshift(Psij)));




    %Psig = Psij ./ abs(Psij) .* sqrt(y_obs);
    Psig = sqrt(data(ii).I).*exp(i*angle(Psij));%%flipud(sqrt(data(ii).I)).*exp(i*angle(Psij));%
    grad = Psij - Psig; %%%% NEW % Psij - Psig;
    grad = fftshift(ifftn(fftshift(grad)));
    grad = repmat(grad, [1 1 size(probe,3)]);
    grad = grad .* conj(Qterm);
    grad = grad .* conj(bmtemp);
    gradtot  = 2*grad; %add factor of 2 to make sure that gradient is equivalent


end
