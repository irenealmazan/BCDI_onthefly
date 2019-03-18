function display_data_gen(xy,qx,qy,dir_file,numb)
%jclark
if exist('numb') == 0,numb='';end

lw=1.5;

font_size=25;

xy=log10(xy);


fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

imagesc(qx,qy,xy)
set(gca,'FontSize',round(0.8*font_size))
xlabel('Spatial Frequency (nm^{-1})','FontSize', font_size), ylabel('Spatial Frequency (nm^{-1})','FontSize', font_size)
print(fh, '-dpng','-r300', [dir_file,'ROI_X-Y',numb]);
exportfig(fh,[dir_file,'ROI_X-Y',numb],'Color','rgb','Renderer','zbuffer')



end

