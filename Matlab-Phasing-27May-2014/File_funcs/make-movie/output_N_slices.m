function [ output_args ] = output_N_slices(array,N_slices,orient,sz,save_dir,save_name,ph_range)
%jclark
%output N 2d slices from a 3d array

%check for directory
if isdir(save_dir) ~= 1,mkdir(save_dir);end

if exist('ph_range') ~= 1,ph_range=[];end
%if exist('hsvim') ~= 1,hsvim = 1;end

%set an roi
if sz ~= 1
    xx=[sz(1),sz(2)];
    yy=[sz(3),sz(4)];
    zz=[sz(5),sz(6)];
else
    sz=size(array);
    xx=[1,sz(2)];
    yy=[1,sz(1)];
    zz=[1,sz(3)];   
end
    
    
switch orient
    
    case 'xy'
        aa=yy;
        bb=xx;
        nc=ceil(size(array,3)/2);
    case 'xz'
        aa=zz;
        bb=xx;
        nc=ceil(size(array,1)/2);
    case 'zy'
        aa=yy;
        bb=zz;
        nc=ceil(size(array,2)/2);
end


for qq=1:N_slices
   
    fh = figure(76) ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    ph = extract_3D_slice(array,orient,nc-round(N_slices/2)+qq);

    if isreal(ph) ~= 1,ph=c2image(ph);end
    
    switch orient
        case 'xy'
            imagesc(ph(yy(1):yy(2),xx(1):xx(2),:));
        case 'xz'
            imagesc(ph(zz(1):zz(2),xx(1):xx(2),:));
        case 'zy'
            imagesc(ph(yy(1):yy(2),zz(1):zz(2),:));  
    end
    
    if ndims(ph) == 2,h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');end
    axis equal
    axis tight
    
    if isempty(ph_range) ~= 1,caxis(ph_range);end
    
    nums = check_numlength(qq,2);
    
    print(fh, '-dpng','-r300', [save_dir,save_name,'-',nums]);
end


end

