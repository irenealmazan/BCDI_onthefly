function output_resolution(params,save_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
lw=1.25;
font_size=15;

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

subplot(2,2,1)
plot(params.xline,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Amplitude','FontSize', font_size)
title('X line-out')
axis([0 max(size(params.xline)) min(params.xline)  1.1*max(params.xline)] )

subplot(2,2,2)
plot(params.dx,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Derivative','FontSize', font_size)
title('Derivative of line-out')
axis([0 max(size(params.dx)) 1.1*min(params.dx)  1.1*max(params.dx)] )

subplot(2,2,3)
plot(params.dx_e1,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Derivative','FontSize', font_size)
title('Derivative of line-out (left edge)')
hold on;
plot(params.fit_dx_e1,'LineWidth',.9*lw,'Color','red')
hold off;
%axis([0 max(size(params.dx_e1)) 1.1*min(params.dx_e1)  1.1*max(params.dx_e1)] )

subplot(2,2,4)
plot(-params.dx_e2,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Derivative','FontSize', font_size)
title('Derivative of line-out (right edge)')
hold on;
plot(-params.fit_dx_e2,'LineWidth',.9*lw,'Color','red')
hold off;
%axis([0 max(size(params.dx_e2)) 1.1*min(-params.dx_e2)  1.1*max(-params.dx_e2)] )
try
    saveas(fh, [save_dir,'Resn-X'], 'epsc'); 
    print(fh, '-dpng','-r300', [save_dir,'Resn-X']);
end
%%
fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

subplot(2,2,1)
plot(params.yline,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Amplitude','FontSize', font_size)
title('Y line-out')
axis([0 max(size(params.yline)) min(params.yline)  1.1*max(params.yline)] )

subplot(2,2,2)
plot(params.dy,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Derivative','FontSize', font_size)
title('Derivative of line-out')
axis([0 max(size(params.dy)) min(params.dy)  1.1*max(params.dy)] )

subplot(2,2,3)
plot(params.dy_e1,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Derivative','FontSize', font_size)
title('Derivative of line-out (left edge)')
hold on;
plot(params.fit_dy_e1,'LineWidth',.9*lw,'Color','red')
hold off;
%axis([0 max(size(params.dy_e1)) 1.1*min(params.dy_e1)  1.1*max(params.dy_e1)] )


subplot(2,2,4)
plot(-params.dy_e2,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Derivative','FontSize', font_size)
title('Derivative of line-out (right edge)')
hold on;
plot(-params.fit_dy_e2,'LineWidth',.9*lw,'Color','red')
hold off;
%axis([0 max(size(params.dy_e2)) 1.1*min(-params.dy_e2)  1.1*max(-params.dy_e2)] )
try
    saveas(fh, [save_dir,'Resn-Y'], 'epsc'); 
    print(fh, '-dpng','-r300', [save_dir,'Resn-Y']);
end

%%
fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

subplot(2,2,1)
plot(params.zline,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Amplitude','FontSize', font_size)
title('Z line-out')
axis([0 max(size(params.zline)) min(params.zline)  1.1*max(params.zline)] )


subplot(2,2,2)
plot(params.dz,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Derivative','FontSize', font_size)
title('Derivative of line-out')
axis([0 max(size(params.dz)) min(params.dz)  1.1*max(params.dz)] )

subplot(2,2,3)
plot(params.dz_e1,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Derivative','FontSize', font_size)
title('Derivative of line-out (left edge)')
hold on;
plot(params.fit_dz_e1,'LineWidth',.9*lw,'Color','red')
hold off;
%axis([0 max(size(params.dz_e1)) 1.1*min(params.dz_e1)  1.1*max(params.dz_e1)] )


subplot(2,2,4)
plot(-params.dx_e2,'LineWidth',lw,'Color','blue')
set(gca,'FontSize',round(0.8*font_size))
xlabel('Pixels','FontSize', font_size), ylabel('Derivative','FontSize', font_size)
title('Derivative of line-out (right edge)')
hold on;
plot(-params.fit_dz_e2,'LineWidth',.9*lw,'Color','red')
hold off;
%axis([0 max(size(params.dz_e2)) 1.1*min(-params.dz_e2)  1.1*max(-params.dz_e2)] )
try
    saveas(fh, [save_dir,'Resn-Z'], 'epsc'); 
    print(fh, '-dpng','-r300', [save_dir,'Resn-Z']);
end

%%
try
    fid=fopen([save_dir,'Resn.txt'],'w');
    fprintf(fid,'X sigma \n');
    fprintf(fid,num2str([params.parms_xe1(2)]));
    fprintf(fid,'\n');
    fprintf(fid,num2str([params.parms_xe2(2)]));
    fprintf(fid,'\n Y sigma \n');
    fprintf(fid,num2str([params.parms_ye1(2)]));
    fprintf(fid,'\n');
    fprintf(fid,num2str([params.parms_ye2(2)]));
    fprintf(fid,'\n Z sigma \n');
    fprintf(fid,num2str([params.parms_ze1(2)]));
    fprintf(fid,'\n');
    fprintf(fid,num2str([params.parms_ze2(2)]));
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    try
        fprintf(fid,'X pixel (nm) \n');
        fprintf(fid,num2str([params.dsx]));
        fprintf(fid,'\n');
        fprintf(fid,'Y pixel (nm) \n');
        fprintf(fid,num2str([params.dsy]));
        fprintf(fid,'\n');
        fprintf(fid,'Z pixel (nm) \n');
        fprintf(fid,num2str([params.dsz]));
    end
  
    fclose(fid);
end
end

