% this script plts the geometry of the different momentum transfer vectors:



figure(3);clf;


gbn=1;

% read the position of the detector position:
 Delta = specscan.motor_positions(1);
 Gamma = specscan.motor_positions(6);
 camdist = specscan.motor_positions(27)/1000; % in meters
 
params.det_px=55e-6;               %detector pixel size (m)
params.lam=.1358;                   %wavelength (nm)
params.delta=Delta;                  %delta (degrees)
params.gam=Gamma;                     %gamma (degrees)
params.arm=camdist;                      %camera distance (m)
rockstep = diff(specscan.var1);
params.dth=double(rockstep(1));                      %angular step size
params.dtilt=0;
 
 
Qlabcenter1 = sind(params.delta)*cosd(params.gam)*(2*pi)/params.lam; %x component, inboard/outboard/horizontal
Qlabcenter2 = sind(params.gam)*(2*pi)/params.lam; %y component, vertical/gravity
Qlabcenter3 = (cosd(params.delta)*cosd(params.gam)-1.0)*(2*pi)/params.lam; %z component, along beam
 
ki=[0,-sind(15),sind(90-15)]*(2*pi)/params.lam;
kf= [sind(params.delta)*cosd(params.gam)*(2*pi)/params.lam, sind(params.gam)*(2*pi)/params.lam, cosd(params.delta)*cosd(params.gam)*(2*pi)/params.lam];
 
ki_o = ki;
kf_o = kf;

quiver3(0,0,0, ki(1), ki(2), ki(3), 'r');
hold on;
quiver3(0,0,0, kf(1), kf(2), kf(3), 'k');

Qlab = [Qlabcenter1,Qlabcenter2,Qlabcenter3]';

quiver3(0,0,0, Qlab(1), Qlab(2), Qlab(3), 'r'); 


Qlabcomp = -ki+kf;


quiver3(0,0,0, Qlabcomp(1), Qlabcomp(2),Qlabcomp(3), 'g');
 
Qlab =Qlabcomp;
 
Qlabcenter1 = Qlabcomp(1);
Qlabcenter2 = Qlabcomp(2);
Qlabcenter3 = Qlabcomp(3);











% kf = [0 0 1];
% ki = [0 0 1];
% kmag = 2*pi/lam;
% 
% Ry = [cosd(-del) 0 sind(-del);
%     0 1 0;
%     -sind(-del) 0 cosd(-del)];
% 
% ki = (Ry * ki.').';
% 
% Rx=[1 0 0;
%      0 cosd(gam) -sind(gam);
%      0 sind(gam) cosd(gam)]; 
% ki = (Rx * ki.').';
% 
% qbragg = kf-ki;
% 
% clf; hold on;
% h=di(probe, -.5, 'g',X,Y,Z); alpha(h, .5);
% h=di(NW, -.5, 'y', X,Y,Z); alpha(h,.5);
% %h=di(support,-.5,'b', X,Y,Z); alpha(h, .2);
% %h=di(tube, -.5, 'c', X,Y,Z); alpha(h, .2);
% 
% % K=convhulln(corners);
% % T=delaunayn(corners,{'Qt','Qbb','Qc','Qz'});
% % p=trisurf(K, corners(:,1), corners(:,2), corners(:,3));
% % set(p,'FaceColor','red','EdgeColor','black');
% % alpha(.3);

% quiver3(0,0,0, ki(1), ki(2), ki(3), 'r');
% quiver3(0,0,0, kf(1), kf(2), kf(3), 'k');
% quiver3(0,0,0, qbragg(1), qbragg(2), qbragg(3), 'b');
% [Xd Yd] = meshgrid([-.1 .1]);
% surf(Xd,Yd,ones(size(Xd)));
% 
% %scatter3(pos_probe_exp(:,1), pos_probe_exp(:,2), pos_probe_exp(:,3), 50, Ga_chan, 'filled');
% 
% hold off
% xlabel('x');ylabel('y'); zlabel('z');
% view(-2,53);
% axis image
% colormap hot
% 
% ki_o = ki;
% kf_o = kf;
% qbragg_o = qbragg;
% 
% %%
% %th scan parameters
% 
% dqlist = zeros(numel(fly2Danglist), 3);
% 
% 
% for ii=1:numel(fly2Danglist)
%     
%     dth = fly2Danglist(ii) - thBragg;
%     
%     
%     Ry = [cosd(-dth) 0 sind(-dth);
%         0 1 0;
%         -sind(-dth) 0 cosd(-dth)];
%     
%     ki = (Ry * ki_o.').';
%     kf = (Ry * kf_o.').';
%     
%     dqlist(ii,:) = (kf-ki)-qbragg;
% 
%     hold on;
%     quiver3(0,0,0, ki(1), ki(2), ki(3), 'r--');
%     quiver3(0,0,0, kf(1), kf(2), kf(3), 'k--');
%     %hold off
% end    
% 
% %%
% %scale all results by kmag:
% ki_o = ki_o*kmag; 
% kf_o = kf_o*kmag; 
% qbragg = qbragg_o*kmag;
% dqlist = dqlist * kmag;
% 
% %%% plot the frame:
% %%{
%     [X_square,Y_square,Z_square] = meshgrid([-Npix/2 Npix/2-1]*d2_bragg,[-Npix/2 Npix/2-1]*d2_bragg,[-depth/2 depth/2-1]*d2_bragg);
% X_square_toplot = X_square(:);
% Y_square_toplot = Y_square(:);
% Z_square_toplot = Z_square(:);
% 
% figure(3);
% hold on;
% scatter3(X_square_toplot,Y_square_toplot,Z_square_toplot)
% %}
