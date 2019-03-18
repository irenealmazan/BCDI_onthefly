function handle = setfigsize( fighandle, xsize, ysize)

pos=get(fighandle, 'position');
if nargin<2
    %xsize=pos(3);
    %ysize=pos(4);
    xsize=210;
    ysize=180;
end

set(fighandle, 'position', [pos(1:2) xsize ysize]);
set(fighandle, 'paperpositionmode', 'auto');
set(gcf,'color','w');
