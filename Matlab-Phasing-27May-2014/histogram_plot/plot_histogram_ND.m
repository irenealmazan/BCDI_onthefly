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

nbins=200;

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

% fh=gcf;
% %% fit gaussian to the histogram
% [gh xx]=fit_gauss_data(h(21:end));
% meann=xx(1)+20;
% sig=xx(2);
% A0=xx(3);
% m=xx(4);
% c=xx(5);
% 
% disp([meann/nbins,sig/nbins])
% 
% amp_th=amp.*shrink_wrap(amp,.2,.05);
% amp_sig=std(amp_th((amp_th > 0)));
% amp_mean=mean(amp_th((amp_th > 0)));
% disp([amp_mean,amp_sig])
% 
% fid=fopen([save_dir,'Hist.txt'],'w');
% fprintf(fid,'mean and std from amp values /n');
% fprintf(fid,num2str([amp_mean,amp_sig]));
% fprintf(fid,'/n mean ans std from gauss-hist fit');
% fprintf(fid,num2str([meann/nbins,sig/nbins]));
% fclose(fid)
% 
% xxx=1:numel(h);   %abisca values
% gg=A0*exp(-0.5*(xxx-meann).^2/sig.^2);%+m*xxx+c;
% 
% h1=gca;
% h2 = axes('Position',get(h1,'Position'));
% plot((0:numel(h)-1)/(numel(h)-1),gg,'LineWidth',3)
% axis([0 1 0 maxy])
% set(h2,'YAxisLocation','Right','Color','none','XTickLabel',[],'FontSize',0.1,'YTickLabel',[])
% %set(h2,'YAxisLocation','left','Color','none','XTickLabel',[],'FontSize',font_size)
% set(h2,'XLim',get(h1,'XLim'),'Layer','top')
% set(gcf,'PaperPositionMode','auto')
% print(fh, '-dpng','-r300', [save_dir,params.prefix,'Hist-fit']);
% exportfig(fh,[save_dir,params.prefix,'Hist-fit'],'Color','rgb','Renderer','zbuffer')
% 
% close(fh)
end