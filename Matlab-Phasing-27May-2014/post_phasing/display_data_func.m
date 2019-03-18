function display_data_func(dir_file,width,rotn,upsamp)
%%jclark
% used to display data and output eps and png files
% this is a 'master copy' with dir_file specifying the reconstruction
% directory to look in.  alternatively, copy display_data.m (this is
% display_data_func.m) and run from the reconstruction directory get it 
% to act 'locally'.
%specify width as a fraction of array size.  use 3 element to specifiy for
%each direction or just one value for all

box_on=0;
one_label=1;

%%
try
    data_dir=dir_file;
catch
    name_of_this_file='display_data_func';
    dir_file=which(name_of_this_file);
    dir_file=dir_file(1:findstr(dir_file,name_of_this_file)-1);
    data_dir=dir_file;
end

try
    width;
    if isempty(width) width=.35;end
catch
    width=0.35;
end

try
    rotn;
    if isempty(rotn) rotn=0;end
catch
    rotn=0;
end

if exist('upsamp') ~= 1,upsamp=[];end

if numel(width) ~= 1
    wx=width(1);
    wy=width(2);
    wz=width(3);
else
    wx=width(1);
    wy=width(1);
    wz=width(1);
end

%%
pfile=rdir([dir_file,'**/*PARAMS.mat']);  %find the params files
temp=pfile(1).name;
load(temp)

counter=0;
qq=0;

while counter == 0      %searches for params files that used the .spe's
    qq=1+qq;
    temp=pfile(qq).name;
    load(temp);
    temp_n=char(params.files(1));
    counter=strcmp(lower(temp_n(end-3:end)),'.tif');
    if qq == max(size(pfile)),counter=1;end
end
load(temp) 

params.return_orig_size=1;
params.do_2D=0;
%%

try
    data_dir=params.data_dir;
end

dir_file=[dir_file,'Data-pics/'];       %where to save output images
if isdir(dir_file) ==0,mkdir(dir_file);end  %create the dir if it doesn't exist

bin=params.binning; %get the binning, want to display with minimal binning though

bin(1)=min([bin(1),2]);
bin(2)=min([bin(2),2]); %set a maximimum binning of 2

try 
    aliens=params.aliens;
catch
    aliens=[];
end

files=params.files;
back=params.back;
min_data=params.min_data;

try
    nnc=params.nnc;             % not support yet, will be used for initial cropping
catch
    nnc=[0];
end

%%
full_files=strcat(data_dir,files);
if numel(back) ~= 0,full_bg=strcat(data_dir,back);else full_bg=[];end

data=bin_crop_center(full_files,full_bg,bin,min_data,aliens,nnc);

data=data/max(max(max(data)));

if rotn ~= 0
    data=rot3darb(data,rotn);
    rstr='-R';
else
    rstr=''; 
end

nn=size(data);
nn=[nn(2),nn(1),nn(3)];

%make them even
new_n=[round(width*nn(1)),round(width*nn(2)),round(width*nn(3))];
new_n=( mod(new_n,2) == 1)+new_n;

data=extract_max(data,new_n(1),new_n(2),new_n(3));

sz=size(data);

cent=round(sz/2);

xc=cent(2);
yc=cent(1);
zc=cent(3);

lam=params.lam;
zd=params.arm;
det_px=params.det_px;
dth=params.dth;

if dth == 0,dth=params.dtilt;end

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

qx=(2*pi/lam)*qx*(det_px/zd);

qy=(2*pi/lam)*qy*(det_py/zd);

qz=(2*pi/lam)*qz*(det_pz/zd);

lw=1.5;

font_size=25;

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(qx,x,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
ylabel('Intensity (Arb. units, log scale)','FontSize', font_size)
axis([min(min([qx,-qx])) max(max([qx,-qx])) min(x) .1 ])
if box_on == 0,box off,end
saveas(fh, [dir_file,'line_out_X',rstr], 'epsc');
print(fh, '-dpng','-r300', [dir_file,'line_out_X',rstr]);


fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(qy,y,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
if one_label~=1, ylabel('Intensity (Arb. units, log scale)','FontSize', font_size);end
axis([min(min([qy,-qy])) max(max([qy,-qy])) min(y) .1])
if box_on == 0,box off,end
saveas(fh, [dir_file,'line_out_Y',rstr], 'epsc'); 
print(fh, '-dpng','-r300', [dir_file,'line_out_Y',rstr]);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(qz,z,'LineWidth',lw,'Color','blue') 
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
if one_label~=1, ylabel('Intensity (Arb. units, log scale)','FontSize', font_size);end
axis([min(min([qz,-qz])) max(max([qz,-qz])) min(z) .1])
if box_on == 0,box off,end
saveas(fh, [dir_file,'line_out_Z',rstr], 'epsc'); 
print(fh, '-dpng','-r300', [dir_file,'line_out_Z',rstr]);

%% 

xy=extract_3D_slice(data,'xy',zc);
xz=extract_3D_slice(data,'xz',yc);
yz=extract_3D_slice(data,'yz',xc);

if isempty(upsamp) ~= 1
   xy=zero_pad_ver3(fftshift(fftn(fftshift(xy))),upsamp*size(xy,2),upsamp*size(xy,1));
   xy=fftshift(ifftn(fftshift(xy)));
   xy=real(xy);
end

xy=log10(xy);
xz=log10(xz);
yz=log10(yz);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
%set(fh,'Position',[100,100,1200,1200])
imagesc(qx,qy,xy)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
print(fh, '-dpng','-r300', [dir_file,'ROI_X-Y',rstr]);
exportfig(fh,[dir_file,'ROI_X-Y',rstr],'Color','rgb','Renderer','zbuffer')
%saveas(fh, [dir_file,'ROI_X'],'psc2');
%set(gcf, 'PaperPositionMode', 'auto');
%set(gcf, 'renderer', 'painters');
%print(fh,'-depsc','-r600',[dir_file,'ROI_X']);
%saveas(fh,[dir_file,'ROI_X'],'-depsc')


fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
imagesc(qx,qz,xz)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
print(fh, '-dpng','-r300', [dir_file,'ROI_X-Z',rstr]);
exportfig(fh,[dir_file,'ROI_X-Z',rstr],'Color','rgb','Renderer','zbuffer')
%saveas(fh, [dir_file,'ROI_Y'],'epsc');
%set(gcf, 'renderer', 'painters');
%set(gcf,'Position',[100,100,1200,1200])
%print(fh,'-depsc','-r600',[dir_file,'ROI_Y']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
imagesc(qy,qz,yz)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
print(fh, '-dpng','-r300', [dir_file,'ROI_Y-Z',rstr]);
exportfig(fh,[dir_file,'ROI_Y-Z',rstr],'Color','rgb','Renderer','zbuffer')
%saveas(fh, [dir_file,'ROI_Z'],'epsc');
%set(gcf, 'renderer', 'painters');
%set(gcf,'Position',[100,100,1200,1200])
%print(fh,'-depsc','-r600',[dir_file,'ROI_Z']);
close all
1;
end

