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

support_edge = facet_spacing * edgepad;

%make parallel planes parallel to each pair of facets
v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
v2 = [corners(1,:) - corners(8,:)]; v2 = v2/norm(v2);
v3 = cross(v1,v2); v3 = v3/norm(v3);
T1 = v3(1)*X + v3(2)*Y + v3(3)*Z;
T1 = (T1>-support_edge/2 & T1<support_edge/2); %two parallel lines

v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
v2 = [corners(3,:) - corners(8,:)]; v2 = v2/norm(v2);
v3 = cross(v1,v2); v3 = v3/norm(v3);
T2 = v3(1)*X + v3(2)*Y + v3(3)*Z;
T2 = (T2>-support_edge/2 & T2<support_edge/2); %two parallel lines

v1 = [corners(4,:) - corners(9,:)]; v1 = v1/norm(v1);
v2 = [corners(3,:) - corners(9,:)]; v2 = v2/norm(v2);
v3 = cross(v1,v2); v3 = v3/norm(v3);
T3 = v3(1)*X + v3(2)*Y + v3(3)*Z;
T3 = (T3>-support_edge/2 & T3<support_edge/2); %two parallel lines

v3 = [0 1 0];
T4 = v3(1)*X + v3(2)*Y + v3(3)*Z;
T4 = (T4>-support_edge/2 & T4<support_edge/2); %two parallel lines


support = T1&T2&T3&T4;

%support = T3;
support = abs(NW);

%sigblur = .007;
%g= 1/(sigblur*sqrt(2*pi)) * exp(-(X.^2 + Y.^2 + Z.^2)/(2*sigblur^2));
%support = fftshift(ifftn( fftn(g).*fftn(NW))) >.1;

rho_guess = rand(Npix,Npix,depth).*0.01 .* exp(i*2*pi*rand(Npix,Npix,depth));
%rho_guess = NW;
%rho_guess = rho;

hold on;
%di(NW,-.5,'g', X,Y,Z);
%h=di(T1, -.5, 'y', X,Y,Z); alpha(h, .5);
h=di(support, -.5, 'y', X,Y,Z); alpha(h, .5);
%plot3(S(1:100:end,1), S(1:100:end,2), S(1:100:end,3), 'k.');
axis image;

K=convhulln(corners);
T=delaunayn(corners,{'Qt','Qbb','Qc','Qz'});
p=trisurf(K, corners(:,1), corners(:,2), corners(:,3));
%p=trimesh(K, corners(:,1), corners(:,2), corners(:,3));
set(p,'FaceColor','red','EdgeColor','black');
%alpha(.3);
%plot3(corners([4 3 9],1), corners([4 3 9],2), corners([4 3 9],3), 'bo');
hold off

xlabel('x'); ylabel('y'); zlabel('z');
drawnow;

%%

Niter_rho = 2000;
Niter_pos = 1;
Niter_theta = 1;
Ncycles = 1;
freq_pos = 10;
freq_rho = 1;
display(['set # rho iterations to ' num2str(Niter_rho) ' temp, position iterations to ' num2str(Niter_pos) ' and cyles ' num2str(Ncycles)])
beta = .8;
rho = rho_guess .* support;
errlist = [];
midsl = round(depth/2);
printind = round( [10:10:100]*(numel(data_exp)/100));
Dmask = (X+Y)>0;
cnt_ntheta = 1;

% prepare the grid and the tools that we need to explore the error as a function of the position of the beam:
ki_norm = ki_o/norm(ki_o);
kf_norm = [0 0 1];

angx = acosd(dot(ki_norm([1 3]), kf_norm([1 3])));
angy = acosd(dot(ki_norm([2 3]), kf_norm([2 3])));

% rotation matrices:

Ry1=[cosd(th) 0 sind(th);
    0 1 0;
    -sind(th) 0 cosd(th)];

Ry2=[cosd(-del) 0 sind(-del);
    0 1 0;
    -sind(-del) 0 cosd(-del)];

Rx=[1 0 0;
    0 cosd(gam) -sind(gam);
    0 sind(gam) cosd(gam)];



% store the initial positions:
 % update data structure & prepare initial guess:
 %{
 for ii = 1:numel(data_exp)
     
     %%% POSITIONS:
     
     %correct positions:     
     %{
     shift_guess(ii,1) = data_exp(ii).xshift_delta;
     shift_guess(ii,2) = data_exp(ii).yshift_delta;
     shift_guess(ii,3) = data_exp(ii).zshift_delta;
     %}
     
     % uncorrect initial positions
     %{
     shift_guess(ii,1) = data_exp(ii).xpos -data_exp(ii).zpos*tand(angx);%data_exp(ii).xshift_delta;
     shift_guess(ii,2) = data_exp(ii).ypos;%data_exp(ii).yshift_delta;
     shift_guess(ii,3) = 0;%data_exp(ii).zshift_delta;
     %}
     
     %{
     dy_guess_ini(ii) = data_exp(ii).ypos;%data_exp(ii).yshift_delta;
     dx_guess_ini(ii) = data_exp(ii).xpos -data_exp(ii).zpos*tand(angx);%data_exp(ii).xshift_delta;
     dz_guess_ini(ii) = 0;%data_exp(ii).zshift_delta;
    %}

    inibeam(ii).A = probe;
    %inibeam(ii).A = circshift(probe, round([shift_guess(ii,2)  shift_guess(ii,1) shift_guess(ii,3)]/d2_bragg));

    %%% THETA POSITIONS
    
    %uncorrect theta positions
    %{
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
  
figure(5); clf; setfigsize(gcf, 1000,500); pause(.1);

for CC = 1:Ncycles
    
      cnt = 1;  % counter for the error surfaces
    
    %PIE ITERATIONS
    %%{
    for nrho = 1:Niter_rho
        
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
                        
            [gPIEiter] = calc_grad_multiangle_Irene(beammat(ii).A, rho, data_exp(ii));
            rho = rho - beta * D* gPIEiter / depth;
            rho = rho .*support;
            if ismember(ii,printind) fprintf('%2.0f%%..', 100*ii/numel(data_exp)); end
%            
           %} 
        end
       
        fprintf('\n');
        
        [err] = calc_error_multiangle_Irene(beammat, rho, data_exp);
        fprintf('     error: %4.4d \n', err);
        errlist = [errlist err];
         %}
        
        subplot(131); imagecomp(rho(:,:,midsl)); colorbar; axis image; %zoom(1.5);
        subplot(132); plot(log10(errlist));
        drawnow;
        toc;
        
        % POSITION ANNEALING
        %{
        tic;
        if mod(nrho,freq_pos) == 0
            
            [delta_guess,err_single,vector_ii_match,err_pos] = NW_posannealing_brutal(Niter_pos,probe,data_exp,rho,Ry1,Ry2,Rx,freq_rho,support,beta,depth);
            
            [err_theta_pos] = NW_annealing_brutal_and_total(Niter_pos,probe,data_exp,rho,Ry1,Ry2,Rx,freq_rho,support,beta_par,depth,dth_disp);
            
            
            shift_guess(:,1) = delta_guess(:,1);
            shift_guess(:,2) = delta_guess(:,2);
            shift_guess(:,3) = delta_guess(:,3);
            
            % update the beam:
             beammat(ii).A = circshift(probe, round([shift_guess(ii,2) shift_guess(ii,1) shift_guess(ii,3)]/d2_bragg));

            errlist = [errlist err_pos];
            number_match_allcycles = [number_match_allcycles sum(vector_ii_match)];
            limit_match = [limit_match numel(data_exp)];
        end
        
        subplot(133); plot(number_match_allcycles); hold on; plot(limit_match);
        drawnow;
        toc;
        %}
        
         % THETA ANNEALING
        %%{
        tic;
        if mod(nrho,freq_pos) == 0
            
         
            [dth_new,grad_final_theta,dq_shift] = NW_theta_annealing_Irene_v3(beammat, rho,data_exp,Niter_theta,index_to_distort);
            %[dth_new,grad_final_theta,dq_shift] = NW_theta_annealing_Irene_v3(inibeam, rho,data_exp, thBragg);
            %[dth_new,dq_shift, fitQerr,err] = NW_theta_annealing_Irene(probe, rho,data_exp,[-0.005:1e-3:0.005], thBragg,shift_guess);
            
            % store the shift
            %theta_guess = dth_new;
            for ii = index_to_distort%1:numel(data_exp)%
               data_exp(ii).dqshift(:) = dq_shift(ii,:); 
               %data_exp(ii).dth_delta_guess(cnt_ntheta) = theta_guess(ii);
               data_exp(ii).dth_new = dth_new(ii);
            end
            
            
            [err] = calc_error_multiangle_Irene(beammat, rho, data_exp);
            %[err] = calc_error_multiangle_Irene(inibeam, rho, data_exp);
            fprintf('     error: %4.4d \n', err);
            errlist = [errlist err];
            figure(5);subplot(132); plot(log10(errlist));
            drawnow;
           % display([' dth guessed = ' num2str(dth_new(index_to_distort)) 'dth true' num2str(data_exp(index_to_distort).dth + data_exp(index_to_distort).dth_delta)])
            
             [mod_ini_square,mod_fin_square] = NW_calculate_final_dq_improvement(data_exp);
            
             mod_fin_square_inter(cnt_ntheta,:) = mod_fin_square;
            cnt_ntheta = cnt_ntheta + 1;
        end
        
        %subplot(133); plot(number_match_allcycles); hold on; plot(limit_match);
        drawnow;
        toc;
        %}
        
        
        
      
       
        
       
    end
    %}
    
   
    
    
    % POSITION ANNEALING ITERATIONS
    %{
    [delta_guess,err_single,vector_ii_match,err_pos] = NW_posannealing_brutal(Niter_pos,probe,data_exp,rho,Ry1,Ry2,Rx);
    
     dx_guess = squeeze(delta_guess(:,1));
     dy_guess = squeeze(delta_guess(:,2));
     dz_guess = squeeze(delta_guess(:,3));
%     
    errlist = [errlist err_pos];
    number_match_allcycles = [number_match_allcycles sum(vector_ii_match)];
    limit_match = [limit_match numel(data_exp)];
    
    
    fprintf('\n');
    
    %[err] = calc_error_multiangle( probe, rho, data_exp);
    fprintf('     error: %4.4d \n', err);
    errlist = [errlist err];
    
    number_match_allcycles = [number_match_allcycles sum(vector_ii_match)];
    limit_match = [limit_match numel(data_exp)];
    
    subplot(131); imagecomp(rho(:,:,midsl)); colorbar; axis image; %zoom(1.5);
    subplot(132); plot(log10(errlist));
    subplot(133); plot(number_match_allcycles); hold on; plot(limit_match);
    drawnow;
    toc;
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



