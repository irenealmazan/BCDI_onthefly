function  sum_align_reconstructions(top_dir,pref)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
try
    top_dir;
catch
    %top_dir='/Users/jesseclark/Documents/MATLAB/data_analysis/Au1110/Moyu/L5150/Crystal2/rand-starts/';
end
%%

norm=1;

recursive0={'CVl','NM','CVgps'};%,
quad=0;     %d them in quadrature
mult=1;
rem_ramp=1;

conj_ref_rec=0; %conj reflect the reference

%%
for qqq=1:numel(recursive0)
    recursive=char(recursive0(qqq));
    save_name=upper(recursive);
    
    try
        pref;
        amp_names={};
        coh_names={};
        ph_names={};
        
        for ww=1:numel(pref),
            amp_f=rdir([top_dir,'*','Rnd','**/*','Rnd',num2str(pref(ww)),'-*',recursive,'*AMP.rec']);
            amp_names=[amp_names,char(amp_f.name)];
            
            ph_f=rdir([top_dir,'*','Rnd','**/*','Rnd',num2str(pref(ww)),'-*',recursive,'*PH.rec']);
            ph_names=[ph_names,char(ph_f.name)];
            
            coh_f=rdir([top_dir,'*','Rnd','**/*','Rnd',num2str(pref(ww)),'-*',recursive,'*COH.rec']);
            coh_names=[coh_names,char(coh_f.name)];
            
        end
        amp_names=char(amp_names);
        ph_names=char(ph_names);
        coh_names=char(coh_names);
    catch
        
        amp_f=rdir([top_dir,'*','Rnd','**/*',recursive,'*AMP.rec']);
        ph_f=rdir([top_dir,'*','Rnd','**/*',recursive,'*PH.rec']);
        coh_f=rdir([top_dir,'*','Rnd','**/*',recursive,'*COH.rec']);
        
        amp_names=char(amp_f.name);
        ph_names=char(ph_f.name);
        coh_names=char(coh_f.name);
    
    end
    
    
    %%
    summed=0;
    sum_ph=0;
    sum_coh=0;
    
    if numel(amp_names) > 0,nfiles=numel(amp_names(:,1));else nfiles=0;end

    for qq=1:nfiles

        
        load(strtrim(amp_names(qq,:)),'-mat')
        amp=array;
        amp=amp/sum(amp(:))*1000;  
        
        load(strtrim(ph_names(qq,:)),'-mat')
        ph=array;
        %ph=ph+pi;   ;%set range between 0-2pi  
        
        
        if size(coh_f,1) ~= 0
            load(strtrim(coh_names(qq,:)),'-mat')
            coh=abs(array);
        else
            coh=0;
        end
        
        %take the ramp off and 'zero' the phase
        if rem_ramp == 1
            F0=ifftshift(fftn(fftshift(amp)));
            F1=ifftshift(fftn(fftshift(amp.*exp(i*ph))));
            [h k l]=register_3d_reconstruction(abs(F0),abs(F1));
            crystal=ifftshift(ifftn(fftshift(sub_pixel_shift(F1,h,k,l))));
            amp=abs(crystal);
            ph=atan2(imag(crystal),real(crystal));
            
            SS=shrink_wrap(amp,.3,.5);  %get just the crystal, i.e very tight support
            avg_ph=mean(ph(find(SS > 0)));
            ph=ph-avg_ph;
        end
        
        if qq == 1,
            big_arr=zeros([size(amp),nfiles]);
            big_ph=big_arr;
            
            if conj_ref_rec == 1
               disp('Conjugating and reflecting the first reconstruction....') 
               cnj = conj_reflect(amp.*exp(i*ph)); 
               amp=abs(cnj);
               ph=atan2(imag(cnj),real(cnj));
            end
            
            ref=amp;    %reference one
            
            big_arr(:,:,:,qq)=amp;    %save them (after checking for cnj)
            big_ph(:,:,:,qq)=ph;
        end
        if qq > 1
            cnj_rf=is_conj_ref(ref,amp);
            
            if cnj_rf ~=0,crystal=conj_reflect(amp.*exp(i*ph));else
                crystal=(amp.*exp(i*ph));end
            
            amp=abs(crystal);
            
            [h k l]=register_3d_reconstruction(ref,amp);

            next_c=(sub_pixel_shift(crystal,h,k,l));
            next=abs(next_c);
            next_ph=atan2(imag(next_c),real(next_c));
            %next=real(sub_pixel_shift(amp,h,k,l));
            %next_ph=real(sub_pixel_shift(ph,h,k,l));

            big_arr(:,:,:,qq)=next;    %save them (after checking for cnj)
        
            big_ph(:,:,:,qq)=next_ph;
        end
        
        if mult == -1
            if qq == 1,sum_coh = 1;end
            sum_coh=sum_coh.*coh;
            if qq == nfiles, sum_coh=sum_coh.^(1.0/nfiles);end
        else
            sum_coh=sum_coh+coh;
        end
            
    end

    if nfiles ~= 0
    
        if quad == 1
            summed=sqrt(sum(big_arr.^2,4));
        elseif mult == 1
            summed=abs(prod(big_arr,4)).^(1.0/nfiles);
        else
            summed=mean(big_arr,4);
        end

        summed_ph=mean(big_ph,4);
        
        save_dir=[top_dir,save_name];
        if isdir(save_dir) == 0,mkdir(save_dir);end

        %% find a params.py file
        php=rdir([top_dir,'**/*phasingparams.py']);
        fid=fopen(char(php(1).name));            %open the file
        str=fscanf(fid,'%c');       %read in the file contents
        fclose(fid); 
        fid=fopen([save_dir,'/','phasingparams.py'],'w');
        fprintf(fid,'%c',str);
        fclose(fid);
        %% find a PARAMS.mat file
        php=rdir([top_dir,'**/*PARAMS.mat']);
        load(char(php(1).name));            %open the file
        save([save_dir,'/','PARAMS.mat'],'params');

        %% output vtk conversion script
        output_python_script(save_dir,'/MATLAB_to_vtk_Qv-AVG.py');

        %% save some png's and eps's
        save_avg_rec(summed,save_dir,'AMP')
        %% calc T form array

        T=ross_tform(params.sz(2),params.sz(1),params.sz(3),params.lam,params.delta,params.gam,...
                     params.det_px*params.binning(1),params.det_px*params.binning(2),params.dth,params.arm); 
        [dsx dsy dsz]=get_sample_pixel_sizes(params.sz(2),params.sz(1),params.sz(3),T);
        %%
        array=summed;

        %array=array/max(abs(array(:)));
        xx=plot_histogram3D(array(find(shrink_wrap(array,.15,.5) > 0) ),[save_dir,'/']);
        %array=array;
        if norm == 1,
            disp('Normalising reconstruction, set norm=0 for no normalisation...')
            array=array/max(array(:));
        end
        save([save_dir,'/',save_name,'-AMP.rec'],'array')
        [ params ] = determine_resolution(array);
        params.dsx=dsx;params.dsy=dsy;params.dsz=dsz;
        output_resolution(params,[save_dir,'/']);

        array=summed_ph;
        save([save_dir,'/',save_name,'-PH.rec'],'array')
        save_avg_rec(array,save_dir,'PH')
        save_avg_rec(array.*shrink_wrap(summed,.3,.3),save_dir,'PH-S')
        
        
        array=std(abs(big_arr),1,4);
        save([save_dir,'/',save_name,'-Std.rec'],'array')
        save_avg_rec(abs(array),save_dir,'Std-AMP')

        avg=summed.*shrink_wrap(summed,.2,.05);
        array=array/mean(avg(avg > 0));
        save([save_dir,'/',save_name,'-Std-mean.rec'],'array')
        save_avg_rec(abs(array),save_dir,'Std-mean-AMP')

        if size(coh_f,1) ~= 0 
            array=sum_coh;
            array=array/max(abs(array(:)));
            save([save_dir,'/',save_name,'-COH.rec'],'array')
        end
    end
    
    disp('Done....')
end

end

function [shifted]=sub_pixel_shift(array,row_shift,col_shift,z_shift)

buf2ft=fftn(array);
[nr,nc,nz]=size(buf2ft);
Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
Nz = ifftshift([-fix(nz/2):ceil(nz/2)-1]);
[Nc,Nr,Nz] = meshgrid(Nc,Nr,Nz);
Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc-z_shift*Nz/nz));
%Greg = Greg*exp(i*diffphase);
shifted=ifftn(Greg);

end

function cnj = conj_reflect(array)

F=ifftshift(fftn(fftshift(array)));

cnj=ifftshift(ifftn(fftshift(conj(F))));


end

function cnj_rf=is_conj_ref(a,b)

c1=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(conj_reflect(b),.1,.1));
c2=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(b,.1,.1));

c1=max(c1(:));
c2=max(c2(:));

if c1 > c2,cnj_rf=1;else cnj_rf=0;end


end

function save_avg_rec(pn,save_dir,prefix)
sx=size(pn);
fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

subplot(2,3,1)
imagesc((pn(:,:,sx(3)/2)))
xlabel('x')
ylabel('y')
colorbar
subplot(2,3,2)
imagesc(squeeze((pn(:,sx(2)/2,:))))
colorbar
xlabel('z')
ylabel('y')
colorbar
subplot(2,3,3)
imagesc(squeeze((pn(sx(1)/2,:,:))))
xlabel('z')
ylabel('x')
colorbar
subplot(2,3,4)
imagesc( squeeze(sum((pn),3)  ))
xlabel('x')
ylabel('y')
colorbar
subplot(2,3,5)
imagesc( squeeze(sum((pn),2)  ))
xlabel('z')
ylabel('y')
colorbar
subplot(2,3,6)
imagesc( squeeze(sum((pn),1)  ))
xlabel('z')
ylabel('x')
ww=0;
colorbar
try
    saveas(fh, [save_dir,'/',prefix], 'epsc');
    print(fh, '-dpng','-r300', [save_dir,'/',prefix]);
end
close(fh);

%%
% SS=shrink_wrap(abs(pn),.3,.5);
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% subplot(2,3,1)
% phase=atan2(imag(pn),real(pn));
% imagesc((phase(:,:,sx(3)/2)))
% xlabel('x')
% ylabel('y')
% colorbar
% subplot(2,3,2)
% imagesc(squeeze((phase(:,sx(2)/2,:))))
% colorbar
% xlabel('z')
% ylabel('y')
% colorbar
% subplot(2,3,3)
% imagesc(squeeze((phase(sx(1)/2,:,:))))
% xlabel('z')
% ylabel('x')
% colorbar
% subplot(2,3,4)
% imagesc((phase(:,:,sx(3)/2)))
% xlabel('x')
% ylabel('y')
% colorbar
% subplot(2,3,5)
% imagesc(squeeze((phase(:,sx(2)/2,:))))
% colorbar
% xlabel('z')
% ylabel('y')
% colorbar
% subplot(2,3,6)
% imagesc(squeeze((phase(sx(1)/2,:,:))))
% xlabel('z')
% ylabel('x')
% colorbar
% try
%     saveas(fh, [save_dir,'/',prefix,'-PH'], 'epsc');
%     print(fh, '-dpng','-r300', [save_dir,'/',prefix,'-PH']);
% end
% close(fh);

end

function  xx=plot_histogram3D(pn,save_dir,maxy )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

amp=abs(pn);

nbins=50;
font_size=25;
lw=1.5;


%amp_th=amp.*shrink_wrap(amp,.2,.05);
[h x]=hist(amp(:),nbins);

try
    maxy;
catch
    maxy=max(h)*1.5;
end


fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
%plot(x,h,'LineWidth',lw,'Color','blue')
bar(x,h,1,'stacked','r');%,'bar_color','b')
set(gca,'FontSize',round(0.8*font_size),'FontName','Helvetica')
xlabel('Amplitude Value','FontSize', font_size,'FontName','Helvetica') 
ylabel('Number of Voxels','FontSize', font_size,'FontName','Helvetica')
axis([0 max(x(:)) 0 maxy])

print(fh, '-dpng','-r300', [save_dir,'Hist']);
exportfig(fh,[save_dir,'Hist'],'Color','rgb','Renderer','zbuffer')


%% fit gaussian to the histogram
[gh xx]=fit_gauss_data(h(21:end));
meann=xx(1)+20;
sig=xx(2);
A0=xx(3);
m=xx(4);
c=xx(5);

disp([meann/nbins,sig/nbins])

amp_th=amp.*shrink_wrap(amp,.2,.05);
amp_sig=std(amp_th((amp_th > 0)));
amp_mean=mean(amp_th((amp_th > 0)));
disp([amp_mean,amp_sig])

fid=fopen([save_dir,'Hist.txt'],'w');
fprintf(fid,'mean and std from amp values /n');
fprintf(fid,num2str([amp_mean,amp_sig]));
fprintf(fid,'/n mean ans std from gauss-hist fit');
fprintf(fid,num2str([meann/nbins,sig/nbins]));
fclose(fid)

xxx=1:numel(h);   %abisca values
gg=A0*exp(-0.5*(xxx-meann).^2/sig.^2);%+m*xxx+c;

h1=gca;
h2 = axes('Position',get(h1,'Position'));
plot((0:numel(h)-1)/(numel(h)-1),gg,'LineWidth',3)
axis([0 1 0 maxy])
set(h2,'YAxisLocation','Right','Color','none','XTickLabel',[],'FontSize',0.1,'YTickLabel',[])
%set(h2,'YAxisLocation','left','Color','none','XTickLabel',[],'FontSize',font_size)
set(h2,'XLim',get(h1,'XLim'),'Layer','top')
set(gcf,'PaperPositionMode','auto')
print(fh, '-dpng','-r300', [save_dir,'Hist-fit']);
exportfig(fh,[save_dir,'Hist-fit'],'Color','rgb','Renderer','zbuffer')

xx(1)=meann/nbins;
xx(2)=sig/nbins;

end

function [gg x] = fit_gauss_data(data)

options = optimset('Display','off','Algorithm','interior-point');
 
x0(1)=mean(data(data > .2*max(data(:))));

x0(2)=std(data(data > .2*max(data(:))));

x0(3)=max(data);

x0(4)=-1/5;
x0(5)=1/20*x0(3);

lb=[0,0,0,-5,-x0(3)];
ub=[numel(data),numel(data),5*x0(3),5,x0(3)];

f=@(x)gauss_fit(x,data);
x=fmincon(f,x0,[],[],[],[],[lb],[ub],[],options);


mn=x(1);
sig=x(2);
A0=x(3);
m=x(4);
c=x(5);

xx=1:numel(data);   %abisca values


gg=A0*exp(-0.5*(xx-mn).^2/sig.^2)+m*xx+c;


end

function E = gauss_fit(x,data,range )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%y=A0*exp(-0.5 x*x/(sig*sig))+mx+c

mean=x(1);
sig=x(2);
A0=x(3);
m=x(4);
c=x(5);


xx=1:numel(data);   %abisca values


gauss=A0*exp(-0.5*(xx-mean).^2/sig.^2)+m*xx+c;

data=data(:);
gauss=gauss(:);

E=sum(abs(gauss-data));


end

function [dsx dsy dsz]=get_sample_pixel_sizes(Nx,Ny,Nz,T)

in=zeros([3,2]);

in(1,:)=(1:2)/Nx;
in(2,:)=(1:2)/Ny;
in(3,:)=(1:2)/Nz;

[ out] = Rot_tform( in,T );

dsx=abs(out(1,1)-out(1,2) );
dsy=abs(out(2,1)-out(2,2) );
dsz=abs(out(3,1)-out(3,2) );

end
