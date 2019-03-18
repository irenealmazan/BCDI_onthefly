function histogram_plot(x,h,maxx,maxy,save_dir,color,xaxis )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
try
    xaxis;
catch
    xaxis='Amplitude Value';
end

font_size=25;
lw=1.5;

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
%plot(x,h,'LineWidth',lw,'Color','blue')
bar(x,h,1,'stacked',color);%,'bar_color','b')
set(gca,'FontSize',round(0.8*font_size),'FontName','Helvetica')

xlabel(xaxis,'FontSize', font_size,'FontName','Helvetica') 

ylabel('Number of Voxels','FontSize', font_size,'FontName','Helvetica')

if numel(maxx) ~= 0 & numel(maxy) ~= 0,axis([0 maxx 0 maxy]);end
box off
print(fh, '-dpng','-r300', [save_dir,'Hist']);
exportfig(fh,[save_dir,'Hist'],'Color','rgb','Renderer','zbuffer')



end

