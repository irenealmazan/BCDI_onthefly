function display_2_coh
%%
dir2='/Volumes/JA-RULE/Users/jesseclark/Documents/MATLAB/data_analysis/Au0710/279_281_283/rand-starts-2k-1/CVL/'
dir1='/Volumes/JA-RULE/Users/jesseclark/Documents/MATLAB/data_analysis/Au0710/Crystal1/261-50-2/rand-starts-2k-1/CVL/'
%dir2='/Users/jesseclark/Documents/MATLAB/data_analysis/Au0710/279_281_283/rand-starts-27x27x27-fftconv-previous/CVL/';
%dir1='/Users/jesseclark/Documents/MATLAB/data_analysis/Au0710/Crystal1/261-50-2/rand-starts/CVL/';


%%
dirs={dir1,dir2};
name_of_this_file='display_2_coh';
dir_file=which(name_of_this_file);
dir_file=dir_file(1:findstr(dir_file,name_of_this_file)-1);
data_dir=dir_file;

colors={'blue','red'};

for ww=1:numel(dirs)
    dir_file=char(dirs(ww));
    pfile=rdir([dir_file,'**/*.csv']);
    sz=size(pfile);
    
    xd=0;
    yd=0;
    zd=0;

    min_c=0.0;  %the minimum to display
    min_g=0.5;   %threhsold for the gauss fit
    lw=1.5;
    font_size=25;
    clr=char(colors(ww));
    norm=1; %normalize
    no_label=1;             %leave label off y and z
    sep_scale=[-1500/2,1500/2];       %do the same scale
    center=1;
    interp=1;               %resample
    plot_fit=1;           %plot the fitted function as dotted line
    
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
                
                qx=x_xx;
                x=x_scalars;
                dx=abs(qx(2)-qx(1));
                
                [gh xx]=fit_gauss_data(x(x > min_g));
                
                lx=sqrt(2*(xx(2)*dx)^2)*sqrt(-log(.88));
                fwhm_x=2.0*sqrt(2*log(2))*xx(2)*dx;
                disp(['lx (.88) = ',num2str(lx),' nm'])
                disp(['lx (HWHM) = ',num2str(fwhm_x/2),' nm'])
                fh = figure(21) ; % returns the handle to the figure object
                if ww ~= 1,hold;end
                set(fh, 'color', 'white'); % sets the color to white 
                
                if center == 1,
                    x=center_array(x);
                    qx=linspace((max(qx)-min(qx))/2,-(max(qx)-min(qx))/2,numel(qx));
                end
                
                plot(qx,x,'LineWidth',lw,'Color',clr)
                set(gca,'FontSize',round(0.8*font_size))
                xlabel('Seperation (nm)','FontSize', font_size) 
                ylabel('Degree of coherence','FontSize', font_size)

                if plot_fit == 1
                   if ww == 1,hold;end
                   plot(qx,exp(-0.5*qx.^2/(dx*xx(2))^2),'--','LineWidth',lw*.8,'Color',clr)
                   hold;
                end
                
                ind=find(x>=min_c);
                if numel(sep_scale) == 2
                    axis([sep_scale(1) sep_scale(2) min_c 1.02 ])
                else
                    axis([min([qx(min(ind))*1.3,-1.3*qx(max(ind))]) max([-1.3*qx(min(ind)),1.3*qx(max(ind))]) min_c 1.02 ])
                end
                if ww == numel(dirs),
                    saveas(fh, [dir_file,'coh_X'], 'epsc');
                    print(fh, '-dpng','-r300', [dir_file,'coh_X']);
                end 
                if ww ~= 1,hold;end
                
            case 'coh_y.csv'
                yd=csvread(dname0,1,0);
                y_scalars=yd(:,1);
                y_xx=yd(:,5);
                if norm == 1,y_scalars=y_scalars/max(y_scalars);end
                
                qy=y_xx;
                y=y_scalars;
                
                if center == 1,
                    y=center_array(y);
                    qy=linspace((max(qy)-min(qy))/2,-(max(qy)-min(qy))/2,numel(qy));
                end
                
                dy=abs(qy(2)-qy(1));
                
                [gh yy]=fit_gauss_data(y(y > min_g));
                ly=sqrt(2*(yy(2)*dy)^2)*sqrt(-log(.88));
                fwhm_y=2.0*sqrt(2*log(2))*yy(2)*dy;
                
                disp(['ly (.88) = ',num2str(ly),' nm'])
                disp(['ly (HWHM) = ',num2str(fwhm_y/2),' nm'])
                fh = figure(31) ; % returns the handle to the figure object
                if ww ~= 1,hold;end
                set(fh, 'color', 'white'); % sets the color to white 
                plot(qy,y,'LineWidth',lw,'Color',clr)
                set(gca,'FontSize',round(0.8*font_size))
                xlabel('Seperation (nm)','FontSize', font_size)
                
                if no_label ~= 1,ylabel('Degree of coherence','FontSize', font_size),else ylabel('Degree of coherence','FontSize', font_size,'Color','white'),end
    
                if plot_fit == 1
                   if ww == 1,hold;end
                   plot(qy,exp(-0.5*qy.^2/(dy*yy(2))^2),'--','LineWidth',lw*.8,'Color',clr)
                   hold;
                end
                
                ind=find(y>=min_c);
                if numel(sep_scale) == 2
                    axis([sep_scale(1) sep_scale(2) min_c 1.02 ])
                else
                    axis([min([qy(min(ind))*1.3,-1.3*qy(max(ind))]) max([-1.3*qy(min(ind)),1.3*qy(max(ind))]) min_c 1.02 ])
                end
                
                if ww == numel(dirs),
                    saveas(fh, [dir_file,'coh_Y'], 'epsc'); 
                    print(fh, '-dpng','-r300', [dir_file,'coh_Y']);
                end
                if ww ~= 1,hold;end
            case 'coh_z.csv'
                zd=csvread(dname0,1,0);
                z_scalars=zd(:,1);
                z_xx=zd(:,6);
                if norm == 1,z_scalars=z_scalars/max(z_scalars);end
                
                qz=z_xx;
                z=z_scalars;
                
                if center == 1,
                    z=center_array(z);
                    qz=linspace((max(qz)-min(qz))/2,-(max(qz)-min(qz))/2,numel(qz));
                end
                dz=abs(qz(2)-qz(1));
                
                [gh zz]=fit_gauss_data(z(z > min_g));
                lz=sqrt(2*(zz(2)*dz)^2)*sqrt(-log(.88));
                fwhm_z=2.0*sqrt(2*log(2))*zz(2)*dz;
                
                disp(['lz (.88) = ',num2str(lz),' nm'])
                
                disp(['lz (HWHM) = ',num2str(fwhm_z/2),' nm'])
                
                fh = figure(41) ; % returns the handle to the figure object
                if ww ~= 1,hold;end
                set(fh, 'color', 'white'); % sets the color to white 
                plot(qz,z,'LineWidth',lw,'Color',clr) 
                set(gca,'FontSize',round(0.8*font_size))
                xlabel('Seperation (nm)','FontSize', font_size) 
                
                if no_label ~= 1,ylabel('Degree of coherence','FontSize', font_size),else ylabel('Degree of coherence','FontSize', font_size,'Color','white'),end
             
                if plot_fit == 1
                   if ww == 1,hold;end
                   plot(qz,exp(-0.5*qz.^2/(dz*zz(2))^2),'--','LineWidth',lw*.8,'Color',clr)
                   hold;
                end
                ind=find(z>=min_c);
                if numel(sep_scale) == 2
                  axis([sep_scale(1) sep_scale(2) min_c 1.02 ])
                else
                  axis([min([qz(min(ind))*1.3,-1.3*qz(max(ind))]) max([-1.3*qz(min(ind)),1.3*qz(max(ind))]) min_c 1.02 ])
                end
                
                if ww == numel(dirs),
                    saveas(fh, [dir_file,'coh_Z'], 'epsc'); 
                    print(fh, '-dpng','-r300', [dir_file,'coh_Z']);
                end
                if ww ~= 1,hold;end
        end

    end
end

% end
end

function [gg x] = fit_gauss_data(data)

data=data(logical(1-isnan(data)));

options = optimset('Display','off','Algorithm','interior-point');
 
x0(1)=0;%mean(data(data > .2*max(data(:))));

ind=find(data > .88);%std(data(data > .2*max(data(:))));
x0(2)=ind(2)-ind(1);

x0(3)=max(data);

%x0(4)=-1/5;
%x0(5)=1/20*x0(3);

lb=[0,0,0];%,-5,-x0(3)];
ub=[numel(data),numel(data),5*x0(3)];%,5,x0(3)];

f=@(x)gauss_fit(x,data);
x=fmincon(f,x0,[],[],[],[],[lb],[ub],[],options);


mn=x(1);
sig=x(2);
A0=x(3);
%m=x(4);
%c=x(5);

xx=1:numel(data);   %abisca values


gg=A0*exp(-0.5*(xx-mn).^2/sig.^2);%+m*xx+c;


end

function E = gauss_fit(x,data,range )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%y=A0*exp(-0.5 x*x/(sig*sig))+mx+c

mn=x(1);
sig=x(2);
A0=x(3);
%m=x(4);
%c=x(5);


xx=1:numel(data);   %abisca values


gauss=A0*exp(-0.5*(xx-mn).^2/sig.^2);%+m*xx+c;

data=data(:);
gauss=gauss(:);

E=sum(abs(gauss-data));


end