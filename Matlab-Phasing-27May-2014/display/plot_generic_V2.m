function fh=plot_generic_V2(x,y,stuff,y2,y3,y4,y5,y6)

AA=blue2red;
AA=circshift(AA,[64,0]);
BB=reverse(AA(16+1:128-16,:));
ind=round(1:(size(BB,1)/(8)):size(BB,1));
ColorSet=BB(ind,:);

try
    y2;
catch
    y2=[];
end
try
    y3;
catch
    y3=[];
end
try
    y4;
catch
    y4=[];
end
try
    y5;
catch
    y5=[];
end

try
    y6;
catch
    y6=[];
end

try
    stuff.legend;
catch
    stuff.legend=[];
end

try
    stuff.legend_loc;
catch
    stuff.legend_loc='Best';
end

try
    lw=stuff.lw;
catch
    lw=1.5;
end
try
    font_size=stuff.font_size;
catch
    font_size=25;
end
try
    xlab=stuff.xlabel;
catch
    xlab='';
end
try
    ylab=stuff.ylabel;
catch
    ylab='';
end
try
    save_name=stuff.save_name;
catch
    save_name=[];
end

try
    stuff.color;
catch
    stuff.color='blue';
end
try
    stuff.color2;
catch
    stuff.color2='red';
end
try
    stuff.color3;
catch
    stuff.color3='green';
end

try
    stuff.color4;
catch
    stuff.color4='orange';
end

try
    stuff.color5;
catch
    stuff.color5='brown';
end

try
    stuff.color6;
catch
    stuff.color6='grey';
end

try
    stuff.n_temps_cold;
catch
    stuff.n_temps_cold=0;
end
%

try
    stuff.dots_only;  %plot data points only
    if numel(stuff.dots_only) == 1,stuff.dots_only=[stuff.dots_only,0];end 
catch
    stuff.dots_only=[0,0];
end
try
    stuff.bold_plot;    %plot in bold
catch
    stuff.bold_plot=0;
end

try
    stuff.error_y;
catch
    stuff.error_y=[];
end

try
    stuff.error_y2;
catch
    stuff.error_y2=[];
end

try
    stuff.error_y3;
catch
    stuff.error_y3=[];
end

try
    stuff.error_y4;
catch
    stuff.error_y4=[];
end

try
    stuff.error_y5;
catch
    stuff.error_y5=[];
end
try
    stuff.error_y6;
catch
    stuff.error_y6=[];
end

try
    stuff.xrange;
catch
    stuff.xrange=[];
end

try
    stuff.yrange;
catch
    stuff.yrange=[];
end

try
    stuff.axis_equal;
catch
    stuff.axis_equal=[];
end

try
    stuff.markersize;
catch
    stuff.markersize=7;
end
markersize=stuff.markersize;
%markersize=7;

try
    stuff.markerfacecolor;
catch
    stuff.markerfacecolor=stuff.color;
end
try
    stuff.markeredgecolor;
catch
    stuff.markeredgecolor=stuff.color;
end

%%%%%%%%%

%
fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white','position',[100,100,650,420]); % sets the color to white 


switch stuff.dots_only(1)
    case 0
        if isempty(stuff.error_y), plot(x,y,'LineWidth',lw,'Color',stuff.color), else errorbar(x,y,stuff.error_y,'LineWidth',lw,'Color',stuff.color);end 
    case 1
        if isempty(stuff.error_y), plot(x,y,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color,'MarkerEdgeColor',stuff.color), else errorbar(x,y,stuff.error_y,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color,'MarkerEdgeColor',stuff.color);end
    case -1
        if isempty(stuff.error_y), plot(x,y,'-o','MarkerSize',markersize,'MarkerFaceColor',stuff.color,'MarkerEdgeColor',stuff.color,'LineWidth',lw,'Color',stuff.color), else errorbar(x,y,stuff.error_y,'-o','MarkerSize',markersize,'MarkerFaceColor',stuff.color,'MarkerEdgeColor',stuff.color,'LineWidth',lw,'Color',stuff.color);end
    
end


set(gca,'FontSize',round(0.8*font_size))
xlabel(xlab,'FontSize', font_size);ylabel(ylab,'FontSize', font_size); 

hold on;
if numel(y2) ~= 0
    switch stuff.dots_only(2)
%         case 0
%             plot(x,y2,'LineWidth',lw,'Color',stuff.color2)
%         case 1
%             plot(x,y2,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color2,'MarkerEdgeColor',stuff.color2)
%         case -1
%             plot(x,y2,'->','MarkerSize',markersize,'MarkerFaceColor',stuff.color2,'MarkerEdgeColor',stuff.color2,'LineWidth',lw,'Color',stuff.color2)
%         
    case 0
        if isempty(stuff.error_y2), plot(x,y2,'LineWidth',lw,'Color',stuff.color2), else errorbar(x,y2,stuff.error_y2,'LineWidth',lw,'Color',stuff.color2);end 
    case 1
        if isempty(stuff.error_y2), plot(x,y2,'s','MarkerSize',markersize,'MarkerFaceColor',stuff.color2,'MarkerEdgeColor',stuff.color2), else errorbar(x,y2,stuff.error_y2,'s','MarkerSize',markersize,'MarkerFaceColor',stuff.color2,'MarkerEdgeColor',stuff.color2);end
    case -1
        if isempty(stuff.error_y2), plot(x,y2,'-s','MarkerSize',markersize,'MarkerFaceColor',stuff.color2,'MarkerEdgeColor',stuff.color2,'LineWidth',lw,'Color',stuff.color2), else errorbar(x,y2,stuff.error_y2,'-s','MarkerSize',markersize,'MarkerFaceColor',stuff.color2,'MarkerEdgeColor',stuff.color2,'LineWidth',lw,'Color',stuff.color2);end
    
    end
end

if numel(y3) ~= 0
    switch stuff.dots_only(3)
%         case 0
%             plot(x,y3,'LineWidth',lw,'Color',stuff.color3)
%         case 1
%             plot(x,y3,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color3,'MarkerEdgeColor',stuff.color3)
%         case -1
%             plot(x,y3,'-<','MarkerSize',markersize,'MarkerFaceColor',stuff.color3,'MarkerEdgeColor',stuff.color3,'LineWidth',lw,'Color',stuff.color3)
%         
    case 0
        if isempty(stuff.error_y3), plot(x,y3,'LineWidth',lw,'Color',stuff.color3), else errorbar(x,y3,stuff.error_y3,'LineWidth',lw,'Color',stuff.color3);end 
    case 1
        if isempty(stuff.error_y3), plot(x,y3,'>','MarkerSize',markersize,'MarkerFaceColor',stuff.color3,'MarkerEdgeColor',stuff.color3), else errorbar(x,y,stuff.error_y3,'>','MarkerSize',markersize,'MarkerFaceColor',stuff.color3,'MarkerEdgeColor',stuff.color3);end
    case -1
        if isempty(stuff.error_y3), plot(x,y3,'->','MarkerSize',markersize,'MarkerFaceColor',stuff.color3,'MarkerEdgeColor',stuff.color3,'LineWidth',lw,'Color',stuff.color3), else errorbar(x,y3,stuff.error_y3,'->','MarkerSize',markersize,'MarkerFaceColor',stuff.color3,'MarkerEdgeColor',stuff.color3,'LineWidth',lw,'Color',stuff.color3);end
    
    
    end
end

if numel(y4) ~= 0
    switch stuff.dots_only(4)
%         case 0
%             plot(x,y4,'LineWidth',lw,'Color',stuff.color4)
%         case 1
%             plot(x,y4,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color4,'MarkerEdgeColor',stuff.color4)
%         case -1
%             plot(x,y4,'-d','MarkerSize',markersize,'MarkerFaceColor',stuff.color4,'MarkerEdgeColor',stuff.color4,'LineWidth',lw,'Color',stuff.color4)
%         
    case 0
        if isempty(stuff.error_y4), plot(x,y4,'LineWidth',lw,'Color',stuff.color4), else errorbar(x,y4,stuff.error_y4,'LineWidth',lw,'Color',stuff.color4);end 
    case 1
        if isempty(stuff.error_y4), plot(x,y4,'<','MarkerSize',markersize,'MarkerFaceColor',stuff.color4,'MarkerEdgeColor',stuff.color4), else errorbar(x,y4,stuff.error_y4,'<','MarkerSize',markersize,'MarkerFaceColor',stuff.color4,'MarkerEdgeColor',stuff.color4);end
    case -1
        if isempty(stuff.error_y4), plot(x,y4,'-<','MarkerSize',markersize,'MarkerFaceColor',stuff.color4,'MarkerEdgeColor',stuff.color4,'LineWidth',lw,'Color',stuff.color4), else errorbar(x,y4,stuff.error_y4,'-<','MarkerSize',markersize,'MarkerFaceColor',stuff.color4,'MarkerEdgeColor',stuff.color4,'LineWidth',lw,'Color',stuff.color4);end
    
    
    end
end

if numel(y5) ~= 0
    switch stuff.dots_only(5)
%         case 0
%             plot(x,y5,'LineWidth',lw,'Color',stuff.color5)
%         case 1
%             plot(x,y5,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color5,'MarkerEdgeColor',stuff.color5)
%         case -1
%             plot(x,y5,'-d','MarkerSize',markersize,'MarkerFaceColor',stuff.color5,'MarkerEdgeColor',stuff.color5,'LineWidth',lw,'Color',stuff.color5)
%         
    case 0
        if isempty(stuff.error_y5), plot(x,y5,'LineWidth',lw,'Color',stuff.color5), else errorbar(x,y5,stuff.error_y5,'LineWidth',lw,'Color',stuff.color5);end 
    case 1
        if isempty(stuff.error_y5), plot(x,y5,'d','MarkerSize',markersize,'MarkerFaceColor',stuff.color5,'MarkerEdgeColor',stuff.color5), else errorbar(x,y5,stuff.error_y5,'d','MarkerSize',markersize,'MarkerFaceColor',stuff.color5,'MarkerEdgeColor',stuff.color5);end
    case -1
        if isempty(stuff.error_y5), plot(x,y5,'-d','MarkerSize',markersize,'MarkerFaceColor',stuff.color5,'MarkerEdgeColor',stuff.color5,'LineWidth',lw,'Color',stuff.color5), else errorbar(x,y5,stuff.error_y5,'-d','MarkerSize',markersize,'MarkerFaceColor',stuff.color5,'MarkerEdgeColor',stuff.color5,'LineWidth',lw,'Color',stuff.color5);end
    
    end
end

if numel(y6) ~= 0
    switch stuff.dots_only(6)
%         case 0
%             plot(x,y6,'LineWidth',lw,'Color',stuff.color6)
%         case 1
%             plot(x,y6,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color6,'MarkerEdgeColor',stuff.color6)
%         case -1
%             plot(x,y6,'-d','MarkerSize',markersize,'MarkerFaceColor',stuff.color6,'MarkerEdgeColor',stuff.color6,'LineWidth',lw,'Color',stuff.color6)
%         
    case 0
        if isempty(stuff.error_y6), plot(x,y6,'LineWidth',lw,'Color',stuff.color6), else errorbar(x,y6,stuff.error_y6,'LineWidth',lw,'Color',stuff.color6);end 
    case 1
        if isempty(stuff.error_y6), plot(x,y6,'d','MarkerSize',markersize,'MarkerFaceColor',stuff.color6,'MarkerEdgeColor',stuff.color6), else errorbar(x,y6,stuff.error_y6,'d','MarkerSize',markersize,'MarkerFaceColor',stuff.color6,'MarkerEdgeColor',stuff.color6);end
    case -1
        if isempty(stuff.error_y6), plot(x,y6,'-d','MarkerSize',markersize,'MarkerFaceColor',stuff.color6,'MarkerEdgeColor',stuff.color6,'LineWidth',lw,'Color',stuff.color6), else errorbar(x,y6,stuff.error_y6,'-d','MarkerSize',markersize,'MarkerFaceColor',stuff.color6,'MarkerEdgeColor',stuff.color6,'LineWidth',lw,'Color',stuff.color6);end
    
    end
end

hold off;

[aa]=axis;
if numel(y2) == 0,y2=y;end
if numel(y3) == 0,y3=y;end
if numel(y4) == 0,y4=y;end
if numel(y5) == 0,y5=y;end
if numel(y6) == 0,y6=y;end

mm=min([min(y6(:)),min(y5(:)),min(y4(:)),min(y3(:)),min(y2(:)),min(y(:))]);
nn=max([max(y6(:)),max(y5(:)),max(y4(:)),max(y3(:)),max(y2(:)),max(y(:))]);

aa(3)=mm-.1*abs(mm-nn);
aa(4)=nn+.1*abs(mm-nn);

if isempty(stuff.xrange) ~= 1,aa(1)=stuff.xrange(1);aa(2)=stuff.xrange(2);end
if isempty(stuff.yrange) ~= 1,aa(3)=stuff.yrange(1);aa(4)=stuff.yrange(2);end

axis([aa(1) aa(2) aa(3) aa(4)])

NumTicks = 5;
L = get(gca,'YTick');
if mod(numel(L),5) == 1,L(2:2:end)=[];end
set(gca,'YTick',L)

if numel(stuff.legend) ~= 0
    legend(stuff.legend,'Location',stuff.legend_loc) %'best'
end
box off
if stuff.bold_plot == 1,set(gca,'FontWeight','bold','LineWidth',.8*lw);end

if isempty(stuff.axis_equal) ~= 1,axis equal;end

if numel(save_name) ~= 0
    %saveas(fh, save_name, 'epsc');
    print(fh, '-dpng','-r600', save_name);
end

end

function string = extract_from_cfit(cfitthing)


names=coeffnames(cfitthing);

vals=coeffvalues(cfitthing);

string='';

for pp=1:size(vals,2)

    string=[string,'[',(names{pp}),' : ',num2str(vals(pp)),']    '];

end


end


