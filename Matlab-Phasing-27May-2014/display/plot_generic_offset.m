function plot_generic_offset(x,y,stuff)

AA=blue2red;
AA=circshift(AA,[64,0]);
BB=reverse(AA(16+1:128-16,:));
ind=round(1:(size(BB,1)/(8)):size(BB,1));
ColorSet=BB(ind,:);


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
    stuff.offset;
catch
    stuff.offset=zeros([1 size(y,2)]);
end

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

map=colormap(fire);

hold on
for qq=1:size(y,2);
    
    xx=x(:,qq);
    yy=y(:,qq)+(qq-1)*stuff.offset(qq);
    stuff.color=map( (qq-1)*round(size(map,1)/size(y,2))+1,:);
    
    switch stuff.dots_only(1)
        case 0
            if isempty(stuff.error_y), plot(xx,yy,'LineWidth',lw,'Color',stuff.color), else errorbar(xx,yy,stuff.error_y(:,qq),'LineWidth',lw,'Color',stuff.color);end 
        case 1
            if isempty(stuff.error_y), plot(xx,yy,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color,'MarkerEdgeColor',stuff.color), else errorbar(xx,yy,stuff.error_y(:,qq),'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color,'MarkerEdgeColor',stuff.color);end
        case -1
            if isempty(stuff.error_y), plot(xx,yy,'-o','MarkerSize',markersize,'MarkerFaceColor',stuff.color,'MarkerEdgeColor',stuff.color,'LineWidth',lw,'Color',stuff.color), else errorbar(xx,yy,stuff.error_y(:,qq),'-o','MarkerSize',markersize,'MarkerFaceColor',stuff.color,'MarkerEdgeColor',stuff.color,'LineWidth',lw,'Color',stuff.color);end

    end

end
hold off    

set(gca,'FontSize',round(0.8*font_size))
xlabel(xlab,'FontSize', font_size);ylabel(ylab,'FontSize', font_size); 



[aa]=axis;
if numel(y2) == 0,y2=y;end
if numel(y3) == 0,y3=y;end
if numel(y4) == 0,y4=y;end
if numel(y5) == 0,y5=y;end

mm=min([min(y3(:)),min(y2(:)),min(y(:))]);
nn=max([max(y3(:)),max(y2(:)),max(y(:))]);

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


