function save_figure_mult(fh,save_dir,fname,seperate_dirs)
%jclark
%save as eps, png and svg given a figure handle

if exist('seperate_dirs') ~= 1,seperate_dirs=1;end

if seperate_dirs == 1
    save_dir_png=[save_dir,'png/'];
    save_dir_eps=[save_dir,'eps/'];
    save_dir_svg=[save_dir,'svg/'];
else
    save_dir_png=[save_dir];
    save_dir_eps=[save_dir];
    save_dir_svg=[save_dir];
end

if isdir(save_dir_png) ~= 1,mkdir(save_dir_png);end
if isdir(save_dir_eps) ~= 1,mkdir(save_dir_eps);end
if isdir(save_dir_svg) ~= 1,mkdir(save_dir_svg);end


saveas(fh, [save_dir_eps,fname],'epsc');
print(fh, '-dpng','-r300', [save_dir_png,fname]);
plot2svg([save_dir_svg,fname,'.svg'],fh)

end

