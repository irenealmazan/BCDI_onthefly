 % NW_ph_retrieval_Irene_v3 tests the combination of PIE with the steepest descent method
% calculating the slope and the step that we need
% to move in order to minize the error metric which drives or PIE
% minimizationalgorithgm
% 2) effect of the number of iterations to perform over the position
% annealing algorithm
%

% construct the reciprocal space grid:

%global d2_bragg

% [Xq,Yq,Zq] = meshgrid([-Npix/2:Npix/2-1]./d2_bragg, ...
%     [-Npix/2:Npix/2-1]./d2_bragg,...
%     [-depth/2:depth/2-1]./d2_bragg);

Npix_x = size(X,1);
Npix_y = size(Y,2);
Npix_z = size(Z,3);


figure(10); clf;


rho_guess = rand(Npix_x,Npix_y,Npix_z).*0.01 .* exp(i*2*pi*rand(Npix_x,Npix_y,Npix_z));

support  = BCDI_make_initial_support(2);

hold on;
%di(NW,-.5,'g', X,Y,Z);
%h=di(T1, -.5, 'y', X,Y,Z); alpha(h, .5);
h=di(support, -.5, 'y', X,Y,Z); alpha(h, .5);
%plot3(S(1:100:end,1), S(1:100:end,2), S(1:100:end,3), 'k.');
axis image;

%%

Niter_rho = 2000;
Niter_pos = 1;
Niter_theta = 1;
Ncycles = 1;
freq_pos = 15;
freq_rho = 1;
freq_support = 20;
freq_HIO = 100;
display(['set # rho iterations to ' num2str(Niter_rho) ' temp, position iterations to ' num2str(Niter_pos) ' and cyles ' num2str(Ncycles)])
beta = .8;
rho = rho_guess .* support;
rho_iter(1).rho = rho;
errlist = [];
midsl = round(depth/2);
printind = round( [10:10:100]*(numel(data_exp)/100));
Dmask = (X+Y)>0;
cnt_ntheta = 1;
ER = repmat([ones(1,15) zeros(1,5)]',100,1);
HIO = repmat([zeros(1,15) ones(1,5)]',100,1);
support_update_array = zeros(1,Niter_rho);
theta_corr_array = zeros(1,Niter_rho);

max_rock_index = find(rock == max(rock));




% store the initial positions:
 % update data structure & prepare initial guess:
 %%{
 for ii = 1:numel(data_exp)
     
    

    inibeam(ii).A = probe;
    %inibeam(ii).A = circshift(probe, round([shift_guess(ii,2)  shift_guess(ii,1) shift_guess(ii,3)]/d2_bragg));

    %%% THETA POSITIONS
    
    %uncorrect theta positions
    %%{
    dth_initials(ii) = data_exp(ii).dth;
    data_exp(ii).dth_new = data_exp(ii).dth; 
    data_exp(ii).dqshift = [data_exp(ii).dqx data_exp(ii).dqy data_exp(ii).dqz]; 
    %}
    
    %correct theta positions
    %{
     data_exp(ii).dqshift =  data_exp(ii).dqshift_delta; 
    %}
    
    
 end
 
 
 
  
 
[err] = calc_error_multiangle_Irene(inibeam, rho, data_exp);
fprintf('initial  error: %4.4d \n', err);
errlist = [errlist err];

%}
 
 % store the errors:
number_match_allcycles = [];
limit_match = [];
  
figure(5); %clf; setfigsize(gcf, 1000,500); pause(.1);
subplot(151); imagecomp(rho(:,:,midsl)); colorbar; axis image; %zoom(1.5);
subplot(152); plot(log10(errlist));
subplot(155);plot(dth_initials,rock,'-xr');title('rocking curve');

drawnow;

for CC = 1:Ncycles
    
      cnt = 1;  % counter for the error surfaces
    
    %PIE ITERATIONS
    %%{
    for nrho = 16:Niter_rho
        
        tic;
        err=0;
       
        %%{
        gPIEiter = 0;
        
        fprintf('PIE iter %i: ', nrho);
        orderrand = randperm(numel(data_exp));
        for ii = 1:numel(data_exp)
            
            beammat(ii).A = probe;
           %beammat(ii).A = circshift(probe, round([shift_guess(ii,2) shift_guess(ii,1) shift_guess(ii,3)]/d2_bragg));
            
           %%{
            D = 1/(max(max(max(abs(beammat(ii).A).^2))));
                        
            
            
%             if mod(nrho,freq_HIO)
%                HIO = 1;
%                ER = 0;
%             else
%                 HIO = 0;
%                 ER = 1;
%             end
            
            if ER(nrho)
                for kk = 1:5
                    [gPIEiter,Psij_calc(ii).Psij] = calc_grad_multiangle_Irene(beammat(ii).A, rho, data_exp(ii));
                    rho = rho - beta * D* gPIEiter / depth;
                end
                %{
                fft_rho = fftshift(fftn(fftshift(rho)));
                Psij_calc(ii).Psij =  sum(fft_rho,3);
                fft_rho2 = sqrt(data_exp(ii).I).*exp(1i*angle(fft_rho));
                rho = fftshift(ifftn(fftshift(fft_rho2)));%rho - beta * D* gPIEiter / depth;
                %}
                
                rho = rho .*support;
            elseif HIO(nrho)
                [gPIEiter,Psij_calc(ii).Psij] = calc_grad_multiangle_Irene(beammat(ii).A, rho, data_exp(ii));
                ind_update = find(support ~= 0);  
                rho(ind_update) = rho(ind_update) - 0.5 * D* gPIEiter(ind_update) / depth;
                ind_update = find(support == 0); 
                 rho(ind_update) = rho(ind_update) - beta*(rho(ind_update) - 0.5 * D* gPIEiter(ind_update) / depth);
            end
            
            % update the support
            if mod(nrho,freq_support) == 0
                [support] = BCDI_update_support(rho);
                rho_iter(nrho).rho = rho;
                support_update_array(nrho) = 1;
            end
            
            if ismember(ii,printind) fprintf('%2.0f%%..', 100*ii/numel(data_exp)); end
%            
            subplot(151); imagecomp(rho(:,:,midsl)); colorbar; axis image; %zoom(1.5);
            %subplot(142); plot(log10(errlist));
            subplot(153); imagesc(sqrt(data_exp(ii).I)); axis image;title(['Experimental dp ' num2str(ii)]);
            subplot(154); imagesc(sqrt(Psij_calc(ii).Psij.*conj(Psij_calc(ii).Psij)));axis image;title('Calculated dp');
            drawnow;
            
           %} 
        end
       
        fprintf('\n');
        
        [err] = calc_error_multiangle_Irene(beammat, rho, data_exp);
        fprintf('     error: %4.4d \n', err);
        errlist = [errlist err];
         %}
        
        subplot(151); imagecomp(rho(:,:,midsl)); colorbar; axis image; %zoom(1.5);
        subplot(152); plot(log10(errlist));
        subplot(153); imagesc(sqrt(data_exp(max_rock_index).I)); axis image;title('Experimental dp');
        subplot(154); imagesc(sqrt(Psij_calc(max_rock_index).Psij.*conj(Psij_calc(max_rock_index).Psij)));axis image;title('Calculated dp');
        drawnow;
        toc;
        
  
         % THETA ANNEALING
        %%{
        tic;
        if mod(nrho,freq_pos) == 0
             threshold = 0.1;
            index_to_distort = find(rock> threshold*rock(max_rock_index));
            [dth_new,grad_final_theta,dq_shift] = NW_theta_annealing_Irene_v3(beammat, rho,data_exp,Niter_theta,index_to_distort',thBragg2/2,cnt_ntheta);
          
            
            % store the shift
            %theta_guess = dth_new;
            for ii = index_to_distort'%1:numel(data_exp)%
               data_exp(ii).dqshift(:) = dq_shift(ii,:); 
               data_exp(ii).dth_new_iter(cnt_ntheta) = dth_new(ii);
               data_exp(ii).dth_new = dth_new(ii);
            end
            
            
            [err] = calc_error_multiangle_Irene(beammat, rho, data_exp);
            %[err] = calc_error_multiangle_Irene(inibeam, rho, data_exp);
            fprintf('     error: %4.4d \n', err);
            errlist = [errlist err];
            figure(5);subplot(152); plot(log10(errlist));
            subplot(155);plot(dth_new,rock,'-xr');title('rocking curve');
            drawnow;
           % display([' dth guessed = ' num2str(dth_new(index_to_distort)) 'dth true' num2str(data_exp(index_to_distort).dth + data_exp(index_to_distort).dth_delta)])
           %if mod(cnt_ntheta,5) ==0 
           
             theta_corr_array(nrho) = 1;
             %[mod_ini_square,mod_fin_square] = NW_calculate_final_dq_improvement(data_exp);
            
             %mod_fin_square_inter(cnt_ntheta,:) = mod_fin_square;
            %cnt_ntheta = cnt_ntheta + 1;
        end
        
        %subplot(133); plot(number_match_allcycles); hold on; plot(limit_match);
        drawnow;
        toc;
        %}
               
       
    end
    %}
    
   
    
    
    
end

errlist = [errlist calc_error_multiangle(probe, rho, data_exp)];

return

%%
figure(6); clf;
%0.05 for cut = 100 (1.4)
for ii = 0.3
    clf;
    rhotemp = rho.*exp(i*ii*Z/d2_bragg);
    imagesc(squeeze(angle(rhotemp(130,:,:))))
    %imagecomp(squeeze(rhotemp(100,:,:)), [],[],[],.3);
    axis image
    colormap hsv
    colorbar
    title(['ii = ' num2str(ii)])
    pause(0.3)
    
end
%%
for kk = -0.1:0.01:0.1
    clf;
    rhotemp2 = rhotemp.*exp(1i*kk);
    imagesc(squeeze(angle(rhotemp2(188,:,:))));
    %imagecomp(squeeze(rhotemp2(100,:,:)), [],[],[],.3);
    axis image
    colormap hsv
    colorbar
    title(['kk = ' num2str(kk)])
    pause(0.3)
end
%%
%0 for cut = 100 (1.4)
for jj = -.4:0.01:-.1
    clf;
    rhotemp2 = rhotemp.*exp(i*jj*X/d2_bragg);
    %imagesc(squeeze(angle(rhotemp2(100,:,:))));
    imagecomp(squeeze(rhotemp2(188,:,:)), [],[],[],.3);
    axis image
    colormap hsv
    colorbar
    title(['jj = ' num2str(jj)])
    pause(0.3)
end

%%

figure(401);clf; setfigsize(gcf, 1500,500);
l = mean(rho(:,[-10:10]+Npix/2),2);
subplot(211);
plot(imax, angle(l));
subplot(212);
plot(imax, abs(l));

%%
% Compare reconstructions
%%{
reconvals = struct;
recons = [0, 1, 2, 3, 4, 5];
addpath(genpath('\..\..\m_scripts'))
addpath(genpath('\Reconstructions'))
for ll = 1:numel(recons)
    load(['Reconstructions\NW_recon_2110_it25_edgepad1' num2str(recons(ll)) '.mat']);
    reconvals(ll).rho = rho;
    reconvals(ll).err = errlist;
end
load('Reconstructions\NW_recon_2110_it96_edgepad11.mat')
reconvals(7).rho = rho;
reconvals(7).err = errlist;

%compare error
figure(404); clf;
leglabel = [];
for ii = 1:numel(recons)
    semilogy(reconvals(ii).err,'LineWidth',2); hold on;
    leglabel = [leglabel; 'edgepad = 1.' num2str(recons(ii))];
end
semilogy(reconvals(7).err,'--','LineWidth',2)
legend(leglabel);

%Compare imagecomp of rho
figure(403); clf; setfigsize(gcf, 1000,500);
axes('position', [0 0 .5 1]);
for ii = 1:depth
    subplot(1,5,1)
    imagecomp(rho1(:,:,ii), [],[],[],.5);
    colormap hsv;
    subplot(1,5,2)
    imagecomp(rho2(:,:,ii), [],[],[],.5);
    colormap hsv;
    subplot(1,5,3)
    imagecomp(rho3(:,:,ii), [],[],[],.5);
    colormap hsv;
    subplot(1,5,4)
    imagecomp(rho4(:,:,ii), [],[],[],.5);
    colormap hsv;
    subplot(1,5,5)
    imagecomp(rho5(:,:,ii), [],[],[],.5);
    colormap hsv;
    pause(0.4);clf;
end


% compare cuts in Amplitude and Phase
figure(401);clf;
for ii = 1:numel(recons)+1
    reconvals(ii).phaseline = angle(mean(reconvals(ii).rho(:,[-10:10]+Npix/2,depth/2),2));
    reconvals(ii).ampline = abs(mean(reconvals(ii).rho(:,[-10:10]+Npix/2,depth/2),2));
end

for ii = 1:(7)
    subplot(7,2,(ii*2-1)); hold on; plot(Y(:,1,1), reconvals(ii).phaseline);
    if ii == 7
        ylabel('1.1, 97 it')
    else
        ylabel(['1.' num2str(recons(ii))])
    end
    if ii == 1
        title('Phase')
    end
    subplot(7,2,ii*2); hold on; plot(Y(:,1,1), reconvals(ii).ampline);
    if ii == 7
        ylabel('1.1, 97 it')
    else
        ylabel(['1.' num2str(recons(ii))])
    end
    if ii == 1
        title('Amplitude')
    end
end



