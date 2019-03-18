function [gradtot, Psij] = calc_grad_multiangle_Irene(probe, rho, data)

global d2_bragg X Y Z

gradtot = zeros(size(rho));

for ii=1:numel(data)

%     Qterm = exp(i* data(ii).dqx * X) .* ...
%             exp(i* data(ii).dqy * Y) .* ...
%             exp(i* data(ii).dqz * Z);

    Qterm = exp(i* data(ii).dqshift(1) * X) .* ...
        exp(i* data(ii).dqshift(2) * Y) .* ...
        exp(i* data(ii).dqshift(3) * Z);
    
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
    grad = 2*grad; %add factor of 2 to make sure that gradient is equivalent

    gradtot = gradtot + grad;
end
