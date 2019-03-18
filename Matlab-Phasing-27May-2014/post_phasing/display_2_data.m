%%
dir1='L5150/Crystal2/';            %directory of the first data
dir2='L5650/Crystal2/';        %directory of the second data
dirs={dir1,dir2};

output_data=0;                  %save data (will still save lineout) 0=no
%%
name_of_this_file='display_2_data';
dir_file=which(name_of_this_file);
dir_file=dir_file(1:findstr(dir_file,name_of_this_file)-1);
data_dir=dir_file;
nnn=max(size(dirs));

for qq=1:nnn
    
    file_path=[dir_file,char(dirs(qq))];
    pfile=rdir([file_path,'**/*PARAMS.mat']);
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
    min_data=0;%params.min_data;

    nnc=[1,1,1];             % not support yet, will be used for initial cropping
    full_files=strcat(file_path,files);
    if numel(back) ~= 0,full_bg=strcat(file_path,back);else full_bg=[];end

    data=bin_crop_center(full_files,full_bg,bin,min_data,aliens,nnc);

    data=data/max(max(max(data)));

    nn=size(data);
    nn=[nn(2),nn(1),nn(3)];
    
    if qq == 1,crn=[round(0.35*nn(1)),round(0.35*nn(2)),round(0.35*nn(3))];end
    data=extract_max(data,crn(1),crn(2),crn(3));

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

    qx=(2*pi/lam)*qx*(det_px/zd);

    qy=(2*pi/lam)*qy*(det_py/zd);

    qz=(2*pi/lam)*qz*(det_pz/zd);

    lw=1.5;

    font_size=25;
    
    if qq == 1                  %create array to store the lineouts
        line_outs_x=zeros([nnn,max(size(x))]);
        line_outs_y=zeros([nnn,max(size(y))]);
        line_outs_z=zeros([nnn,max(size(z))]);
    end    
    
    line_outs_x(qq,:)=x;
    line_outs_y(qq,:)=y;
    line_outs_z(qq,:)=z;
    
        
    if output_data == 1
        
        xy=extract_3D_slice(data,'xy',zc);
        xz=extract_3D_slice(data,'xz',yc);
        yz=extract_3D_slice(data,'yz',xc);

        xy=log10(xy);
        xz=log10(xz);
        yz=log10(yz);
    
        fh = figure ; % returns the handle to the figure object
        set(fh, 'color', 'white'); % sets the color to white 
        %set(fh,'Position',[100,100,1200,1200])
        imagesc(qx,qy,xy)
        set(gca,'FontSize',round(0.8*font_size))
        xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
        print(fh, '-dpng','-r300', [dir_file,'ROI_X',num2str(qq)]);
        exportfig(fh,[dir_file,'ROI_X',num2str(qq)],'Color','rgb','Renderer','zbuffer')
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
        print(fh, '-dpng','-r300', [dir_file,'ROI_Y',num2str(qq)]);
        exportfig(fh,[dir_file,'ROI_Y',num2str(qq)],'Color','rgb','Renderer','zbuffer')
        %saveas(fh, [dir_file,'ROI_Y'],'epsc');
        %set(gcf, 'renderer', 'painters');
        %set(gcf,'Position',[100,100,1200,1200])
        %print(fh,'-depsc','-r600',[dir_file,'ROI_Y']);

        fh = figure ; % returns the handle to the figure object
        set(fh, 'color', 'white'); % sets the color to white 
        imagesc(qy,qz,yz)
        set(gca,'FontSize',round(0.8*font_size))
        xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
        print(fh, '-dpng','-r300', [dir_file,'ROI_Z',num2str(qq)]);
        exportfig(fh,[dir_file,'ROI_Z',num2str(qq)],'Color','rgb','Renderer','zbuffer')
        %saveas(fh, [dir_file,'ROI_Z'],'epsc');
        %set(gcf, 'renderer', 'painters');
        %set(gcf,'Position',[100,100,1200,1200])
        %print(fh,'-depsc','-r600',[dir_file,'ROI_Z']);

    end
    %% 

  
end



fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(qx,line_outs_x(1,:),'LineWidth',lw,'Color','blue');hold
plot(qx,line_outs_x(2,:)-2,'LineWidth',lw,'Color','blue');hold

set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Intensity (Arb. units, log scale)','FontSize', font_size)
axis([min(min([qx,-qx])) max(max([qx,-qx])) min(x-2) .1 ])
saveas(fh, [dir_file,'line_out_X-2'], 'epsc');
print(fh, '-dpng','-r300', [dir_file,'line_out_X-2']);


fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(qy,line_outs_y(1,:),'LineWidth',lw,'Color','blue');hold
plot(qy,line_outs_y(2,:)-2,'LineWidth',lw,'Color','blue');hold
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Intensity (Arb. units, log scale)','FontSize', font_size,'color','white')
axis([min(min([qy,-qy])) max(max([qy,-qy])) min(y-2) .1])
saveas(fh, [dir_file,'line_out_Y-2'], 'epsc'); 
print(fh, '-dpng','-r300', [dir_file,'line_out_Y-2']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
plot(qz,line_outs_z(1,:),'LineWidth',lw,'Color','blue');hold
plot(qz,line_outs_z(2,:)-2,'LineWidth',lw,'Color','blue');hold
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Intensity (Arb. units, log scale)','FontSize', font_size,'color','white')
axis([min(min([qz,-qz])) max(max([qz,-qz])) min(z-2) .1])
saveas(fh, [dir_file,'line_out_Z-2'], 'epsc'); 
print(fh, '-dpng','-r300', [dir_file,'line_out_Z-2']);