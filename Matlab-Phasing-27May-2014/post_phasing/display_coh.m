%%


%%
name_of_this_file='display_coh';
dir_file=which(name_of_this_file);
dir_file=dir_file(1:findstr(dir_file,name_of_this_file)-1);
data_dir=dir_file;


pfile=rdir([dir_file,'**/*.csv']);
sz=size(pfile);

xd=0;
yd=0;
zd=0;

min_c=0.7;  %the minimum to display

norm=1; %normalize
no_label=1;             %leave label off y and z
sep_scale=[-400,400];       %do the same scale

for qq=1:sz(1)
    
    dname0=pfile(qq).name;
    dname=dname0(end-8:end);
    
    disp(dname0)
    
    switch dname
        
        case 'coh_x.csv'
            xd=csvread(dname0,1,0);
            x_scalars=xd(:,1);
            x_xx=xd(:,4);
            if norm == 1,x_scalars=x_scalars/max(x_scalars);end
            %x_scalars(x_scalars < min_c)=[];
            %x_xx(x_scalars < min_c)=[];
        case 'coh_y.csv'
            yd=csvread(dname0,1,0);
            y_scalars=yd(:,1);
            y_xx=yd(:,5);
            if norm == 1,y_scalars=y_scalars/max(y_scalars);end
            %y_scalars(y_scalars < min_c)=[];
            %y_xx(y_scalars < min_c)=[];
        case 'coh_z.csv'
            zd=csvread(dname0,1,0);
            z_scalars=zd(:,1);
            z_xx=zd(:,6);
            if norm == 1,z_scalars=z_scalars/max(z_scalars);end
            %z_scalars(z_scalars < min_c)=[];
            %z_xx(z_scalars < min_c)=[];
    end
    
end

%%
%function a=  display_3d_data(data,params )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

lw=1.5;
font_size=25;

if numel(xd) > 1,

    qx=x_xx;
    x=x_scalars;
    
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    plot(qx,x,'LineWidth',lw,'Color','blue')
    set(gca,'FontSize',round(0.8*font_size))
    xlabel('Seperation (nm)','FontSize', font_size), ylabel('Degree of coherence','FontSize', font_size)
    
    %axis([min(min([qx,-qx])) max(max([qx,-qx])) min_c 1.02 ])
    ind=find(x>=min_c);
    if numel(sep_scale) == 2
        axis([sep_scale(1) sep_scale(2) min_c 1.02 ])
    else
        axis([min([qx(min(ind))*1.3,-1.3*qx(max(ind))]) max([-1.3*qx(min(ind)),1.3*qx(max(ind))]) min_c 1.02 ])
    end
    %title('x')
    saveas(fh, [dir_file,'coh_X'], 'epsc'); 
    print(fh, '-dpng','-r300', [dir_file,'coh_X']);
end

if numel(yd) > 1,

    qy=y_xx;
    y=y_scalars;
    
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    plot(qy,y,'LineWidth',lw,'Color','blue')
    set(gca,'FontSize',round(0.8*font_size))
    xlabel('Seperation (nm)','FontSize', font_size) 
    
    if no_label ~= 1,ylabel('Degree of coherence','FontSize', font_size),else ylabel('Degree of coherence','FontSize', font_size,'Color','white'),end
    
    %axis([min(min([qy,-qy])) max(max([qy,-qy])) min(y) 1.02])
    ind=find(y>=min_c);
    
    if numel(sep_scale) == 2
        axis([sep_scale(1) sep_scale(2) min_c 1.02 ])
    else
        axis([min([qy(min(ind))*1.3,-1.3*qy(max(ind))]) max([-1.3*qy(min(ind)),1.3*qy(max(ind))]) min_c 1.02 ])
    end
    %title('y')
    saveas(fh, [dir_file,'coh_Y'], 'epsc'); 
    print(fh, '-dpng','-r300', [dir_file,'coh_Y']);
end


if numel(zd) > 1,

    qz=z_xx;
    z=z_scalars;
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    plot(qz,z,'LineWidth',lw,'Color','blue') 
    set(gca,'FontSize',round(0.8*font_size))
    xlabel('Seperation (nm)','FontSize', font_size), 
    
   
    if no_label ~= 1,ylabel('Degree of coherence','FontSize', font_size),else ylabel('Degree of coherence','FontSize', font_size,'Color','white'),end
    ind=find(z>=min_c);
    
    if numel(sep_scale) == 2
        axis([sep_scale(1) sep_scale(2) min_c 1.02 ])
    else
        axis([min([qz(min(ind))*1.3,-1.3*qz(max(ind))]) max([-1.3*qz(min(ind)),1.3*qz(max(ind))]) min_c 1.02 ])
    end
    %title('z')
    saveas(fh, [dir_file,'coh_Z'], 'epsc'); 
    print(fh, '-dpng','-r300', [dir_file,'coh_Z']);
end

