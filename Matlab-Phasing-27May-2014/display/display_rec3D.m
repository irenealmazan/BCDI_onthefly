function  display_rec3D(pn,support,chi,error_DM,save_dir )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sx=size(pn);

%[ array xyz] = center_array( support );
%pn=circshift(real(pn),xyz)+1i*circshift(imag(pn),xyz);
%support=circshift(support,xyz);
%array=0;

[pn yxz]=center_array(pn);
support=circshift(support,yxz);

xyz=center_of_mass(abs(pn).*support);

if ndims(support) == 3
    xyz=-1*round([xyz(2),xyz(1),xyz(3)]);
else
    xyz=-1*round([xyz(2),xyz(1)]);    
end

pn=circshift(real(pn),xyz)+1i*circshift(imag(pn),xyz);
support=circshift(support,xyz);
%

%xyz=center_of_mass(abs(pn).*support);

%if ndims(pn) == 3,xyz=-1*round([xyz(2),xyz(1),xyz(3)]);end
%if ndims(pn) == 2,xyz=-1*round([xyz(2),xyz(1)]);end

%pn=circshift(pn,xyz);%circshift(real(pn),xyz)+1i*circshift(imag(pn),xyz);
%%
%figure

if ndims(pn) == 3
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    subplot(2,4,1)
    p1=loglog((1:numel(chi)),(chi));

    set(p1,'Color','blue','LineWidth',1.0)
    xlabel('iteration')
    ylabel('\chi')
    title('\chi vs iteration number')
    %set(get(gcf,'CurrentAxes'),'FontName','Ariel','FontSize',14)
    %set(get(gcf,'CurrentAxes'), 'FontWeight', 'bold');

    subplot(2,4,5)

    if numel(error_DM) >= 5,p2=loglog( 5:numel(error_DM), error_DM(5:numel(error_DM)));
    else p2=loglog(1:numel(error_DM),error_DM);end

    set(p2,'Color','blue','LineWidth',1.0)
    xlabel('iteration')
    ylabel('||pn+1 - pn||')
    title('iterate error')

    subplot(2,4,2)
    imagesc(abs(pn(:,:,sx(3)/2)))

    subplot(2,4,3)
    imagesc(squeeze(abs(pn(:,sx(2)/2,:))))

    subplot(2,4,4)
    imagesc(squeeze(abs(pn(sx(1)/2,:,:))))

    %subplot(2,4,5)
    %imagesc(squeeze(abs(pn(sx(1)/2,:,:))))

    subplot(2,4,6)
    imagesc( squeeze(sum(abs(pn),3)  ))

    subplot(2,4,7)
    imagesc( squeeze(sum(abs(pn),2)  ))

    subplot(2,4,8)
    imagesc( squeeze(sum(abs(pn),1)  ))

    ww=0;
    
    try
        saveas(fh, [save_dir,'Rec'], 'epsc');
        print(fh, '-dpng','-r300', [save_dir,'Rec']);
    end

end



if ndims(pn) == 2
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    set(fh,'OuterPosition',[200,200,900,600]);
    subplot(2,2,1)
    p1=loglog((1:numel(chi)),(chi));

    set(p1,'Color','blue','LineWidth',1.0)
    xlabel('iteration')
    ylabel('\chi')
    title('\chi vs iteration number')
    %set(get(gcf,'CurrentAxes'),'FontName','Ariel','FontSize',14)
    %set(get(gcf,'CurrentAxes'), 'FontWeight', 'bold');

    subplot(2,2,3)

    if numel(error_DM) >= 5,p2=loglog( 5:numel(error_DM), error_DM(5:numel(error_DM)));
    else p2=loglog(1:numel(error_DM),error_DM);end

    set(p2,'Color','blue','LineWidth',1.0)
    xlabel('iteration')
    ylabel('||pn+1 - pn||')
    title('iterate error')

    subplot(2,2,2)
    imagesc(abs(pn))

    subplot(2,2,4)
    imagesc(abs(support))


    ww=0;
    
end



    
end

