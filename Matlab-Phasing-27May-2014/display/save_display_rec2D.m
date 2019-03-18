function  save_display_rec2D(pn,support,save_dir )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sx=size(pn);


xyz=center_of_mass(abs(pn).*support);
xyz=-1*round([xyz(2),xyz(1)]);
pn=circshift(pn,xyz);%circshift(real(pn),xyz)+1i*circshift(imag(pn),xyz);

%figure
fh = figure(34) ; % returns the handle to the figure object

%set(fh,'OuterPosition',[200,200,1200,600]);
set(fh, 'color', 'white'); % sets the color to white 

%subplot(1,3,1)
imagesc(abs(pn))
xlabel('x')
ylabel('y')

try
    saveas(fh, [save_dir,'Amp'], 'epsc');
    print(fh, '-dpng','-r300', [save_dir,'Amp']);
end
close(fh)  

%%
fh = figure(35) ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

ph=atan2(imag(pn),real(pn));

%subplot(1,3,2)
imagesc((ph))

xlabel('x')
ylabel('y')
try
    saveas(fh, [save_dir,'Ph'], 'epsc');
    print(fh, '-dpng','-r300', [save_dir,'Ph']);
end
close(fh)

%%

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

%subplot(1,3,3)
imagesc((support))

xlabel('x')
ylabel('y')

try
    saveas(fh, [save_dir,'Sup'], 'epsc');
    print(fh, '-dpng','-r300', [save_dir,'Sup']);
end

%try
%    saveas(fh, [save_dir,'Amp+Ph'], 'epsc');
%    print(fh, '-dpng','-r300', [save_dir,'Amp+Ph']);
%end
close(fh)  



end


