
figure(3);clf;


kf = [0 0 1];
ki = [0 0 1];
kmag = 2*pi/lam;

Ry = [cosd(-del) 0 sind(-del);
    0 1 0;
    -sind(-del) 0 cosd(-del)];

ki = (Ry * ki.').';

Rx=[1 0 0;
     0 cosd(gam) -sind(gam);
     0 sind(gam) cosd(gam)]; 
ki = (Rx * ki.').';

qbragg = kf-ki;

clf; hold on;
%h=di(probe, -.5, 'g',X,Y,Z); alpha(h, .5);

%h=di(support,-.5,'b', X,Y,Z); alpha(h, .2);
%h=di(tube, -.5, 'c', X,Y,Z); alpha(h, .2);

% K=convhulln(corners);
% T=delaunayn(corners,{'Qt','Qbb','Qc','Qz'});
% p=trisurf(K, corners(:,1), corners(:,2), corners(:,3));
% set(p,'FaceColor','red','EdgeColor','black');
% alpha(.3);

quiver3(0,0,0, ki(1), ki(2), ki(3), 'r');
quiver3(0,0,0, kf(1), kf(2), kf(3), 'k');
quiver3(0,0,0, qbragg(1), qbragg(2), qbragg(3), 'b');
[Xd Yd] = meshgrid([-.1 .1]);
surf(Xd,Yd,ones(size(Xd)));

%scatter3(pos_probe_exp(:,1), pos_probe_exp(:,2), pos_probe_exp(:,3), 50, Ga_chan, 'filled');

hold off
xlabel('x');ylabel('y'); zlabel('z');
view(-2,53);
axis image
colormap hot

ki_o = ki;
kf_o = kf;
qbragg_o = qbragg;

%%
%th scan parameters

dqlist = zeros(numel(thscanlist), 3);

for ii=1:numel(thscanlist)
    
    dth = thscanvals(ii);
    
    Ry = [cosd(-dth) 0 sind(-dth);
        0 1 0;
        -sind(-dth) 0 cosd(-dth)];
    
    ki = (Ry * ki_o.').';
    kf = (Ry * kf_o.').';
    
    dqlist(ii,:) = (kf-ki)-qbragg;

 %   hold on;
%    quiver3(0,0,0, ki(1), ki(2), ki(3), 'r--');
  %  quiver3(0,0,0, kf(1), kf(2), kf(3), 'k--');
   % hold off
end    

%%
%scale all results by kmag:
ki_o = ki_o*kmag; 
kf_o = kf_o*kmag; 
qbragg = qbragg_o*kmag;
dqlist = dqlist * kmag;

% do the meshgrid:

th_pixel_size = 2*pi/abs((dqlist(end,3) - dqlist(1,3)));

[X Y Z] = meshgrid([-Npix/2:Npix/2-1]*d2_bragg, ...
                    [-Npix/2:Npix/2-1]*d2_bragg,...
                    [-depth/2:depth/2-1]*th_pixel_size);


 % do the object:
 NW_make_InGaAs_nocoreshell_BCDI;
                
%%% plot the frame:
%%{
[X_square,Y_square,Z_square] = meshgrid([-Npix/2 Npix/2-1]*d2_bragg,[-Npix/2 Npix/2-1]*d2_bragg,[-depth/2 depth/2-1]*d2_bragg);
X_square_toplot = X_square(:);
Y_square_toplot = Y_square(:);
Z_square_toplot = Z_square(:);

figure(3);
hold on;
scatter3(X_square_toplot,Y_square_toplot,Z_square_toplot)
h=di(NW, -.5, 'y', X,Y,Z); alpha(h,.5);
%}
