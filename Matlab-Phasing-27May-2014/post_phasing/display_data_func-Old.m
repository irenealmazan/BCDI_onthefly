function display_data_func(dir_file)
%%
% used to display data and output eps and png files
% this is a 'master copy' with dir_file specifying the reconstruction
% directory to look in.  alternatively, copy display_data.m (this is
% display_data_func.m) and run from the reconstruction directory get it 
% to act 'locally'.

%%
try
    dir_file;
catch
    name_of_this_file='display_data_func';
    dir_file=which(name_of_this_file);
    dir_file=dir_file(1:findstr(dir_file,name_of_this_file)-1);
    data_dir=dir_file;
end

pfile=rdir([dir_file,'**/*PARAMS.mat']);
pfile=pfile(1).name;
load(pfile)

%%
bin=[2,2];

try 
    aliens=params.aliens;
catch
    aliens=[];
end

files=params.files;
back=params.back;
min_data=params.min_data;

nnc=[1,1,1];             % not support yet, will be used for initial cropping
full_files=strcat(data_dir,files);
if numel(back) ~= 0,full_bg=strcat(data_dir,back);else full_bg=[];end

data=bin_crop_center(full_files,full_bg,bin,min_data,aliens,nnc);

data=data/max(max(max(data)));

nn=size(data);
nn=[nn(2),nn(1),nn(3)];

data=extract_max(data,round(0.35*nn(1)),round(0.35*nn(2)),round(0.35*nn(3)));



%%
%function a=  display_3d_data(data,params )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

name_of_this_file='display_data';

dir_file=which(name_of_this_file);
dir_file=dir_file(1:findstr(dir_file,name_of_this_file)-1);

sz=size(data);

cent=round(sz/2);

xc=cent(2);
yc=cent(1);
zc=cent(3);

lam=params.lam;
zd=params.arm;
det_px=params.det_px;
dth=params.dth;

try 
    bin;
catch
    bin=params.binning;
end

det_py=det_px*bin(2);
det_px=det_px*bin(1);
det_pz=2*zd*sin(pi/180*dth/2);


x  = extract_1D_slice(data,'x',yc,zc );
x=x/max(abs(x));

y  = extract_1D_slice(data,'y',xc,zc );
y=y/max(abs(y));

z  = extract_1D_slice(data,'z',xc,yc );
z=z/max(abs(z));

x=log10(x);
y=log10(y);
z=log10(z);


maxx=find(x == max(x));

maxy=find(y == max(y));

maxz=find(z == max(z));

if numel(maxx) == 1
    qx=(1:size(x))-maxx;
else 
    qx=(1:size(x))-maxx(1)-0.5;
end

if numel(maxy) == 1
    qy=(1:size(y))-maxy;
else 
    qy=(1:size(y))-maxy(1)-0.5;
end

if numel(maxz) == 1
    qz=(1:size(z))-maxz;
else 
    qz=(1:size(z))-maxz(1)-0.5;
end

hx=lam*zd/det_px./abs(qx);
hx(find(hx == Inf))=zd*lam/det_px;

hy=lam*zd/det_py./abs(qy);
hy(find(hy == Inf))=zd*lam/det_py;

hz=lam*zd/det_pz./abs(qz);
hz(find(hz == Inf))=zd*lam/det_pz;

qx=(4*pi/lam)*qx*(det_px/zd);

qy=(4*pi/lam)*qy*(det_py/zd);

qz=(4*pi/lam)*qz*(det_pz/zd);

lw=1.5;

font_size=25;

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(qx,x,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Intensity (Arb. units, log scale)','FontSize', font_size)
axis([min(min([qx,-qx])) max(max([qx,-qx])) min(x) .1 ])
saveas(fh, [dir_file,'line_out_X'], 'epsc'); 


fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(qy,y,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Intensity (Arb. units, log scale)','FontSize', font_size)
axis([min(min([qy,-qy])) max(max([qy,-qy])) min(y) .1])
saveas(fh, [dir_file,'line_out_Y'], 'epsc'); 

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(qz,z,'LineWidth',lw,'Color','blue') 
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Intensity (Arb. units, log scale)','FontSize', font_size)
axis([min(min([qz,-qz])) max(max([qz,-qz])) min(z) .1])
saveas(fh, [dir_file,'line_out_Z'], 'epsc'); 

%% 

xy=extract_3D_slice(data,'xy',zc);
xz=extract_3D_slice(data,'xz',yc);
yz=extract_3D_slice(data,'yz',xc);

xy=log10(xy);
xz=log10(xz);
yz=log10(yz);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
imagesc(qx,qy,xy)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
print(fh, '-dpng','-r600', [dir_file,'ROI_X']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
imagesc(qx,qz,xz)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
print(fh, '-dpng','-r600', [dir_file,'ROI_Y']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
imagesc(qy,qz,yz)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
print(fh, '-dpng','-r600', [dir_file,'ROI_Z']);

1;
end
