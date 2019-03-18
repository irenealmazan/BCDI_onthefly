function display_2_data_rec()
%%
dir1='Rec-Rnd3-261-ERHIO200-CVl-SW/';            %directory of the first data
dir2='Rec-Rnd3-261-ERHIO200-NM-SW/';
%%
name_of_this_file='display_2_data_rec';
dir_file=which(name_of_this_file);
dir_file=dir_file(1:findstr(dir_file,name_of_this_file)-1);
data_dir=dir_file;


lw=1.5;
font_size=25;



%% load the data and params
file_path=[dir_file,dir1];
[data params]=load_params_data(file_path);


%% load the first one
file_path=[dir_file,dir1];
[pnm_cvl support]=load_reconstruction(file_path);

pns_cvl=pnm_cvl.*support;

data_cvl=abs(fftshift(fftn(ifftshift(pns_cvl))));

[coh]=load_coherence(file_path);

data_cvl_a=sqrt(convnfft(data_cvl.^2,coh,'same'));

%% load the second one
file_path=[dir_file,dir2];
[pnm_nm support]=load_reconstruction(file_path);

pns_nm=pnm_nm.*support;

data_nm=abs(fftshift(fftn(ifftshift(pns_nm))));

%% some prelimanrys
nn=size(data);
nn=[nn(2),nn(1),nn(3)];

crn=[round(0.35*nn(1)),round(0.35*nn(2)),round(0.35*nn(3))];

data=extract_max(data.^2,crn(1),crn(2),crn(3));

[x_m y_m z_m qx_m qy_m qz_m]=get_lineouts_plotting(center_array(data),params,1);

data_nm=extract_max(data_nm.^2,crn(1),crn(2),crn(3));

[x_n y_n z_n qx_n qy_n qz_n]=get_lineouts_plotting(center_array(data_nm),params,1);

data_cvl=extract_max(data_cvl.^2,crn(1),crn(2),crn(3));

[x_c y_c z_c qx_c qy_c qz_c]=get_lineouts_plotting(center_array(data_cvl),params,1);

data_cvl_a=extract_max(data_cvl_a.^2,crn(1),crn(2),crn(3));

[x_ca y_ca z_ca qx_ca qy_ca qz_ca]=get_lineouts_plotting(center_array(data_cvl_a),params,1);

%% plot on a log scale
logs=1;
n_plots=4;
AA=blue2red;
AA=circshift(AA,[64,0]);
BB=reverse(AA(32+1:128-32,:));
ind=round(1:(size(BB,1)/(n_plots)):size(BB,1));
ColorSet=BB(ind,:);

save_name=[dir_file,'x-rec'];
legends={'I^a','I^m','I^u','I^v'};
%legends={'|\psi_{k}^{pc}|^2 \otimes \gamma','I_m','|\psi_{k}^{c}|^2','|\psi_{k}^{pc}|^2'};

plot_line_out(x_n-0,x_m-2,x_ca-3.6,x_c-4.8,qx_m,ColorSet,logs,save_name,legends)

save_name=[dir_file,'y-rec'];
plot_line_out(circshift(y_n-0,1),y_m-2,y_ca-3.6,y_c-4.8,qy_m,ColorSet,0*logs,save_name,[])

save_name=[dir_file,'z-rec'];
plot_line_out(z_n-0,z_m-2,z_ca-3.6,z_c-4.8,qz_m,ColorSet,0*logs,save_name,[])


%% plot on a linear scale
[x_m y_m z_m qx_m qy_m qz_m]=get_lineouts_plotting(data,params,0);
[x_n y_n z_n qx_n qy_n qz_n]=get_lineouts_plotting(data_nm,params,0);
[x_c y_c z_c qx_c qy_c qz_c]=get_lineouts_plotting(data_cvl,params,0);
[x_ca y_ca z_ca qx_ca qy_ca qz_ca]=get_lineouts_plotting(data_cvl_a,params,0);



1;

end

function [pnm support]=load_reconstruction(file_path)

afile=rdir([file_path,'**/*AMP.rec']);
pfile=rdir([file_path,'**/*PH.rec']);
load(afile.name,'-mat')
amp=array;
load(pfile.name,'-mat')
ph=array;

sfile=rdir([file_path,'**/*SUP.rec']);
load(sfile.name,'-mat')
support=array;

pnm=amp.*exp(i*ph);

end

function [data params]=load_params_data(file_path)


pmfile=rdir([file_path,'**/*PARAMS.mat']);
load(pmfile.name,'-mat')

dfile=rdir([file_path,'**/*DATA.rec']);
load(dfile.name,'-mat')
data=array;
array=0;

end

function [coh]=load_coherence(file_path)

pfile=rdir([file_path,'**/*COH.rec']);
load(pfile.name,'-mat')

coh=real(center_array(fftshift(fftn(ifftshift(array)))));


end

function [x y z qx qy qz]=get_lineouts_plotting(data,params,log)

try
    log;
catch
    log=1;
end
    
bin=params.binning;
lam=params.lam;
zd=params.arm;
det_px=params.det_px;
dth=params.dth;
det_py=det_px*bin(2);
det_px=det_px*bin(1);
det_pz=2*zd*sin(pi/180*dth/2);


sz=size(data);

cent=round(sz/2);

xc=cent(2);
yc=cent(1);
zc=cent(3);

x  = extract_1D_slice(data,'x',yc,zc );
x=x/max(abs(x));

y  = extract_1D_slice(data,'y',xc,zc );
y=y/max(abs(y));

z  = extract_1D_slice(data,'z',xc,yc );
z=z/max(abs(z));

if log == 1
    x=log10(x);
    y=log10(y);
    z=log10(z);
end

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


end

function plot_line_out(data1,data2,data3,data4,x,ColorSet,logs,save_name,legends)

lw=1.5;
font_size=25;

mind=min([min(data1(:)),min(data2(:)),min(data3(:)),min(data4(:))]);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(x,data1,'LineWidth',lw,'Color',ColorSet(1,:));hold on
plot(x,data2,'LineWidth',lw,'Color',ColorSet(2,:));
plot(x,data3,'LineWidth',lw,'Color',ColorSet(3,:));
plot(x,data4,'LineWidth',lw,'Color',ColorSet(4,:));
hold off

set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size)

if logs == 1
    ylabel('Intensity (Arb. units, log scale)','FontSize', font_size)
elseif logs == 0
    ylabel('Intensity (Arb. units)','FontSize', font_size,'Color','w')
end
axis([min(min([x,-x])) max(max([x,-x])) mind .1 ])


if numel(legends) ~= 0,legend(legends,'Location','Best');end
if numel(save_name) ~= 0
    saveas(fh, save_name, 'epsc');
    print(fh, '-dpng','-r300', save_name);
end

end


