function [ output_args ] = display_auto_func(data,dir_file,params)
%jclark
%outputs the autocorrelation, assumes its magnitude that is n

lw=1.5;
font_size=25;
pow=.25;

sz=size(data);

try
    lam=params.lam;zd=params.arm;det_px=params.det_px;dth=params.dth;
    if dth ==0,dth=params.dtilt;end
    dx=zd*lam/sz(2)/det_px;dy=zd*lam/sz(1)/det_px;dz=lam/tan(deg2rad(dth*sz(3)));
catch
    dx=1;dy=1;dz=1; 
end

auto=abs(fftshift(ifftn(ifftshift(data.^2)))).^pow;

cent=round(sz/2);
xc=cent(2);
yc=cent(1);
zc=cent(3);
xy=extract_3D_slice(auto,'xy',zc);
xz=extract_3D_slice(auto,'xz',yc);
yz=extract_3D_slice(auto,'yz',xc);

qx=((1:sz(2))-xc)*dx;
qy=((1:sz(1))-yc)*dy;
qz=((1:sz(3))-zc)*dz;

if isdir(dir_file) ~=1,mkdir(dir_file);end

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
imagesc(qx,qy,xy)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Distance (nm)','FontSize', font_size), ylabel('Distance (nm)','FontSize', font_size)
save_figure_mult(fh,dir_file,'ROI_X-Y')

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
imagesc(qx,qz,xz)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Distance (nm)','FontSize', font_size), ylabel('Distance (nm)','FontSize', font_size)
save_figure_mult(fh,dir_file,'ROI_X-Z')

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
imagesc(qy,qz,yz)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Distance (nm)','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
save_figure_mult(fh,dir_file,'ROI_Y-Z')

end

