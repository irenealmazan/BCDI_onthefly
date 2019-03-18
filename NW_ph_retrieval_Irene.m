% construct the reciprocal space grid:

global d2_bragg

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

support = T3&T2&T1;

%support = T3;
%support = abs(NW);

%sigblur = .007;
%g= 1/(sigblur*sqrt(2*pi)) * exp(-(X.^2 + Y.^2 + Z.^2)/(2*sigblur^2));
%support = fftshift(ifftn( fftn(g).*fftn(NW))) >.1;

rho_guess = rand(Npix,Npix,depth)*1e-6 .* exp(i*2*pi*rand(Npix,Npix,depth));
%rho_guess = NW;

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

Niter_rho = 5;
Niter_pos = 1;
Ncycles = 100;
display(['set # rho iterations to ' num2str(Niter_rho) ' temp, position iterations to ' num2str(Niter_pos) ' and cyles ' num2str(Ncycles)])
beta = .8;
rho = rho_guess .* support;
errlist = [];
midsl = round(depth/2);
printind = round( [10:10:100]*(numel(data_exp)/100));
Dmask = (X+Y)>0;

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

%%%% NEW: generate random displacements of the beam with respect to its ideal
% position: 
dx_grid = [-10:5:10]*d2_bragg;
dy_grid = [0]*d2_bragg;
dz_grid = [0]*d2_bragg;

yshift = zeros(numel(data_exp),size(dx_grid,2));
xshift = zeros(numel(data_exp),size(dx_grid,2));
%shift_delta_store = zeros(size(dx_grid,2),3);
 %index = [1 2 2];

% store the initial positions:
 % update data structure & prepare initial guess:
 for ii = 1:numel(data_exp)
     data_exp(ii).ypos_ini = data_exp(ii).ypos;
     data_exp(ii).xpos_ini = data_exp(ii).xpos;
     data_exp(ii).zpos_ini = data_exp(ii).zpos;
     
     dy_guess(ii) = 0;
     dx_guess(ii) = 0;
     dz_guess(ii) = 0;
     
     % calculate initial beam position and error:
     yshift_ini = data_exp(ii).ypos; %This is not generally true, only for horizontal scattering
     xshift_ini = data_exp(ii).xpos - (data_exp(ii).zpos) * tand(angx);
     
     tempbm_ini(ii).A = circshift(probe, round([yshift_ini  xshift_ini 0]/d2_bragg));
 end
 
 
 [err] = calc_error_multiangle_xpos_sim_final(tempbm_ini, rho, data_exp);
 fprintf('    initial  error: %4.4d \n', err);
    
 
 % store the errors:
 err_single = zeros(Ncycles,Niter_pos,numel(data_exp),size(dx_grid,2));
 vector_ii_match = zeros(numel(data_exp),1);
number_match_allcycles = [];
limit_match = [];
  
figure(5); clf; setfigsize(gcf, 1000,500); pause(.1);

for CC = 1:Ncycles
    
    %PIE ITERATIONS
    %%{
    for nrho = 1:Niter_rho
        
        tic;
        err=0;
        
        gPIEiter = 0;
        
        fprintf('PIE iter %i: ', nrho);
        orderrand = randperm(numel(data_exp));
        for ii = 1:numel(data_exp)
            
            if CC == 1
                beammat(ii).A = circshift(probe, round([data_exp(ii).ypos data_exp(ii).xpos data_exp(ii).zpos]/d2_bragg));
            else
                beammat(ii).A = correctbm(ii).A;
            end
            
            D = 1/(max(max(max(abs(beammat(ii).A).^2))));
            
            
            [gPIEiter] = calc_grad_multiangle_Irene(beammat(ii).A, rho, data_exp(ii));
            rho = rho - beta * D* gPIEiter / depth;
            rho = rho .*support;
            if ismember(ii,printind) fprintf('%2.0f%%..', 100*ii/numel(data_exp)); end
            
        end
        fprintf('\n');
        
        [err] = calc_error_multiangle_xpos_sim_final(beammat, rho, data_exp);
        fprintf('     error: %4.4d \n', err);
        errlist = [errlist err];
        
        subplot(131); imagecomp(rho(:,:,midsl)); colorbar; axis image; %zoom(1.5);
        subplot(132); plot(log10(errlist));
        drawnow;
        toc;
        
        
        
    end
    %}
    
   
    % POSITION ANNEALING ITERATIONS
    %%{
    for npos = 1:Niter_pos
        
        for ii = 1:numel(data_exp)
             
            for mm = 1:numel(dx_grid)
                                
                xshift_delta = dx_grid(mm);                
                yshift_delta = 0;
                zshift_delta = 0;
                
                % build the shift vector, to lie in the scaned plane:
                shift_delta = [xshift_delta yshift_delta zshift_delta];
                
                shift_delta = (Ry1*shift_delta')';
                shift_delta = (Ry2*shift_delta')';
                shift_delta = (Rx*shift_delta')';
                
                shift_delta_store(mm,:) = shift_delta;
                
                
                yshift(ii,mm) = data_exp(ii).ypos + shift_delta(2); %This is not generally true, only for horizontal scattering
                xshift(ii,mm) = data_exp(ii).xpos + shift_delta(1) - (data_exp(ii).zpos + shift_delta(3)) * tand(angx);
                
                tempbm = circshift(probe, round([yshift(ii,mm)  xshift(ii,mm) 0]/d2_bragg));
               
                
                [err_single(CC,npos,ii,mm)] = calc_error_multiangle_xpos_sim(tempbm, rho, data_exp(ii));

                
%                 [grad_x(CC,npos,jj,mm)] = calc_grad_proj3d_x(probe,rho,data_exp(ii), [yshift(mm) xshift(mm) 0],Npix_x,Npix_y,Npix_z);
%                
%                 
%                 if mm>1
%                     test_grad_manually;
%                 else
%                     display_grad_x;
%                 end
                
            end
              
            % calculate the true shift to compare:
            yshift_true(ii) = data_exp(ii).ypos + data_exp(ii).yshift_delta; %This is not generally true, only for horizontal scattering
            xshift_true(ii) = data_exp(ii).xpos + data_exp(ii).xshift_delta - (data_exp(ii).zpos + data_exp(ii).zshift_delta) * tand(angx);            
%             
%             tempbm_true = circshift(probe, round([yshift_true(ii)  xshift_true(ii) 0]/d2_bragg));
% 
%             [err_true_position] = calc_error_multiangle_xpos_sim(tempbm_true, rho, data_exp(ii));
%             
             % identifies the index mm which contains the smaller error and stores the displacement in delta_xy
             [minerror,mm_error] = min(err_single(CC,npos,ii,:));
             %delta_yx(CC,jj,:) = [dy_guess(mm_error) dx_guess(mm_error) ];
               
             display(['position ' num2str(ii) ' shift in x = ' num2str(xshift(ii,mm_error)) ', real displacement x = ' num2str(xshift_true(ii)) ])

             % count the number of matches:
              if xshift(ii,mm_error) == xshift_true(ii)                 
                    vector_ii_match(ii) = 1;
              end
             
             % update the beam:
             correctbm(ii).A = circshift(probe, round([yshift(ii,mm_error)  xshift(ii,mm_error) 0]/d2_bragg));

            
            
            %[grad_x(CC,npos,jj)] = calc_grad_proj3d_x(probe,rho, data_exp(ii),[0 dx_guess(ii) 0] ,   Xq,Yq);
            
             %dx_guess(ii) = dx_guess(ii) - grad_x(CC,npos,jj); % see Fineup paper, appendix
            
             % update data structure:
             data_exp(ii).dy_guess = 0;
             data_exp(ii).dx_guess = dx_guess(ii);
             data_exp(ii).dz_guess = 0;
             
             % update data structure:
             data_exp(ii).ypos = data_exp(ii).ypos + data_exp(ii).dy_guess;
             data_exp(ii).xpos = data_exp(ii).xpos + data_exp(ii).dx_guess;
             data_exp(ii).zpos = data_exp(ii).zpos + data_exp(ii).dz_guess;
        end
        
       
       
        fprintf('\n');
        
        [err] = calc_error_multiangle_xpos_sim_final(correctbm, rho, data_exp);
        fprintf('     error: %4.4d \n', err);
        errlist = [errlist err];
        
        number_match_allcycles = [number_match_allcycles sum(vector_ii_match)];
        limit_match = [limit_match numel(data_exp)];
        
        subplot(131); imagecomp(rho(:,:,midsl)); colorbar; axis image; %zoom(1.5);
        subplot(132); plot(log10(errlist));
        subplot(133); plot(number_match_allcycles); hold on; plot(limit_match);
        drawnow;
        toc;
        
        
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



