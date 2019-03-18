function [ prtf ] = Compute_PRTF( rcnst_abs, meas_abs)
%pass in slices of the reconstructed fourier amplitudes and slices of the
%measured fourier amplitudes

%make sure to form rcnst_abs by taking FT(amp>0.2 * exp[i phase*amp>0.2))


[Nx Ny] = size(rcnst_abs); %should only be 2D

%#ok<*BDSCI>
%--------------------------------------------------------------------------
%keep this for the time being, compare speeds for defining it here vs passing it
[xx,yy] = meshgrid(1:Nx,1:Ny);
%--------------------------------------------------------------------------
%the q = 0 value:
jj = 1; kk = 1;

o_diam_circ = (((Nx/2+1-xx)./jj).^2 + ((Ny/2+1-yy)./jj).^2 < 1);

prtf(kk) = sum(sum(o_diam_circ .* rcnst_abs)) / ...
                      (1E-9 + sum(sum(o_diam_circ .* meas_abs)));
%--------------------------------------------------------------------------
jj_max = round(Nx/2); 
annuls_smplng = 1;

%after q = 0, which corresponds to jj = 1, we now need to loop from:
loop_range = ((1 + annuls_smplng) : annuls_smplng : jj_max);
%--------------------------------------------------------------------------
prtf(2 : (1 + length(loop_range))) = 0; 
%--------------------------------------------------------------------------
tic
kk = 2;
for jj = loop_range
    
  o_diam_circ = (((Nx/2+1-xx)./jj).^2 + ((Ny/2+1-yy)./jj).^2 < 1);
  i_diam_circ = (((Nx/2+1-xx)./(jj-1)).^2 + ((Ny/2+1-yy)./(jj-1)).^2 < 1);

  temp22 = o_diam_circ - i_diam_circ;

  %figure(858); imagesc(log10(1+abs( temp22 .* rcnst_abs )))
  
  prtf(kk) = sum(sum(temp22 .* rcnst_abs)) / (1E-9 + sum(sum(temp22 .* meas_abs)));
  %prtf(kk) = sum(sum(temp22 .* meas_abs)) / (1E-9 + sum(sum(temp22 .* rcnst_abs)));
  
  kk = kk+1;

end
toc
%--------------------------------------------------------------------------
figure; semilogx(prtf,'-or','LineWidth',2,'MarkerSize',4); 
xlim([1 length(prtf)])
%ylim([1E-7 1E5])
%--------------------------------------------------------------------------
clear('xx','yy','annuls_smplng','jj_max','loop_range','jj','kk')
clear('temp*','o_diam_circ','i_diam_circ')
%--------------------------------------------------------------------------

end

