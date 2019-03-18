function [ params ] = density_shells_histogram(amp,support,shells,save_dir,params )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

try
    params.maxy;
catch
    params.maxy=1000;
end
try
    params.maxx;
catch
    params.maxx=1;
end
try
    params.xaxis;
catch
    params.xaxis='Amplitude Value';
end

try
    colors=params.colors;
catch
    colors=['b','r','g','m','c','b','r','g','m','c'];
end

try
    params.output_shell_to_csv;
catch
    params.output_shell_to_csv=0;
end
%colors=['b','r','g','m','c'];


%%%%%

shells=abs(get_density_shells(amp,support,shells));
params.shells=shells;      %set the out shell equal to the calc one

hist_dir=[save_dir,'histograms/'];
if isdir(hist_dir) == 0,mkdir(hist_dir);end

for qq=1:size(shells,4)
    
    params.save_dir=[hist_dir,num2str(qq)];
    
    %amp=(amp/max(amp(:)));  %means it will be between 0-1
    
    
    params.color=colors(qq);
    
    plot_histogram_ND(amp((shells(:,:,:,qq) > 0)),params)
    params.shell_vals{qq}=amp((shells(:,:,:,qq) > 0));  %keep values
    
    %output each shell to csv
    if params.output_shell_to_csv==1
       disp('Saving shell values to csv....')
       csvwrite([params.save_dir,params.prefix,'Hist.csv'],amp((shells(:,:,:,qq) > 0)));
    end
    
    if size(shells,4) > 2
        params.save_dir=[hist_dir,'<',num2str(qq)];
        plot_histogram_ND(amp(( sum(shells(:,:,:,qq:size(shells,4)),4) > 0)),params)
    end
    
end


end

function shell_arr=get_density_shells(amp,support,shells)

if sum(support) == 0,support=shrink_wrap(amp,.1,.5);end

inc=1/shells;

levels=linspace(1,0,shells);


amp=amp/max(amp(:));

shell_arr=zeros([size(amp),shells-1]);

for qq=1:shells-1
    
   shell=shrink_support(support,levels(qq));
   
   if qq ~= 1
       shell_arr(:,:,:,qq)=shell;
       shell_arr(:,:,:,qq-1)=shell_arr(:,:,:,qq-1)-shell;
   else
       shell_arr(:,:,:,qq)=shell;
   end
   shell=0;
  
    
end


end

function  plot_histogram_ND(pn,params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
amp=pn;
try
    params.total;
catch
    params.total=[];
end
total=params.total;

if numel(total) > 0,amp=abs(pn)/sum(abs(pn(:))).*total;end

% try
%     params.maxy;
% catch
%     maxy=0.1*numel(find(amp > .2 *max(amp(:))));
%     disp(['maxy - [',num2str(maxy),']'])
%     params.maxy=maxy;
% end
%maxy=params.maxy;

try
    params.maxx;
catch
    params.maxx=max(amp(:));
end
maxx=params.maxx;

save_dir=params.save_dir;

try
    params.color
catch
    params.color='r';
    %color=params.color;
end
try
    params.prefix;
catch
    params.prefix='';
end

nbins=50;

[h x]=imhist(amp(:),nbins);


try
    params.maxy;
catch
    params.maxy=1.1*max(h(:));
end
maxy=params.maxy;

%if numel(params.prefix) == 0,
%    histogram_plot(x,h,maxx,maxy,save_dir,params.color );
%else
histogram_plot(x,h,maxx,maxy,[save_dir,params.prefix],params.color,params.xaxis );
%end

fh=gcf;
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

fid=fopen([save_dir,params.prefix,'Hist.txt'],'w');
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
print(fh, '-dpng','-r300', [save_dir,params.prefix,'Hist-fit']);
exportfig(fh,[save_dir,params.prefix,'Hist-fit'],'Color','rgb','Renderer','zbuffer')

close(fh)
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


