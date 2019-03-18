function  save_display_rec3D(pn,support,save_dir )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sx=size(pn);


xyz=center_of_mass(center_array(abs(pn).*support));

if ndims(pn) == 3
    xyz=-1*round([xyz(2),xyz(1),xyz(3)]);
else
    xyz=-1*round([xyz(2),xyz(1)]);
end

pn=circshift(pn,xyz);%circshift(real(pn),xyz)+1i*circshift(imag(pn),xyz);

%figure
fh = figure(34) ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

if ndims(pn) == 3
    subplot(2,3,1)
    imagesc(abs(pn(:,:,round(sx(3)/2))))
    xlabel('x')
    ylabel('y')

    subplot(2,3,2)
    imagesc(squeeze(abs(pn(:,sx(2)/2,:))))

    xlabel('z')
    ylabel('y')

    subplot(2,3,3)
    imagesc(squeeze(abs(pn(sx(1)/2,:,:))))
    xlabel('z')
    ylabel('x')
    %subplot(2,4,5)
    %imagesc(squeeze(abs(pn(sx(1)/2,:,:))))

    subplot(2,3,4)
    imagesc( squeeze(sum(abs(pn),3)  ))
    xlabel('x')
    ylabel('y')

    subplot(2,3,5)
    imagesc( squeeze(sum(abs(pn),2)  ))
    xlabel('z')
    ylabel('y')

    subplot(2,3,6)
    imagesc( squeeze(sum(abs(pn),1)  ))
    xlabel('z')
    ylabel('x')
    ww=0;
else
    imagesc(abs(pn))
    xlabel('x')
    ylabel('y')
end


try
    saveas(fh, [save_dir,'Amp'], 'epsc');
    print(fh, '-dpng','-r300', [save_dir,'Amp']);
end
close(fh)    

fh = figure(35) ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

ph=atan2(imag(pn),real(pn));

if ndims(pn) == 3

    subplot(2,3,1)
    imagesc((ph(:,:,ceil(sx(3)/2))))

    xlabel('x')
    ylabel('y')

    subplot(2,3,2)
    imagesc(squeeze((ph(:,ceil(sx(2)/2),:))))

    xlabel('z')
    ylabel('y')

    subplot(2,3,3)
    imagesc(squeeze((ph(ceil(sx(1)/2),:,:))))
    xlabel('z')
    ylabel('x')
    %colorbar('location','southoutside')

    subplot(2,3,4)
    imagesc((support(:,:,ceil(sx(3)/2))))

    xlabel('x')
    ylabel('y')

    subplot(2,3,5)
    imagesc(squeeze((support(:,ceil(sx(2)/2),:))))

    xlabel('z')
    ylabel('y')

    subplot(2,3,6)
    imagesc(squeeze((support(ceil(sx(1)/2),:,:))))
    xlabel('z')
    ylabel('x')
    ww=0;
    
else
    imagesc((ph))
    xlabel('x')
    ylabel('y')
end

try
    saveas(fh, [save_dir,'Ph+S'], 'epsc');
    print(fh, '-dpng','-r300', [save_dir,'Ph+S']);
end
close(fh)  



end


