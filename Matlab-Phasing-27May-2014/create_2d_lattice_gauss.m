function [lattice] = create_2d_lattice_gauss(spacing,dsize,dx,dy)
%jclark
%creates a 2d or 3d latttice
%dsize is x size, ysize

%gauss 'atoms'
ax=5;
ay=5;
sx=0.75;
sy=0.75;

npointsx=round(dsize(1)/spacing);
npointsy=round(dsize(2)/spacing);

xp=1:spacing:dsize(1);
yp=1:spacing:dsize(2);

lattice=zeros([dsize(2),dsize(1)]);

for nx=1:numel(xp)
    for ny=1:numel(yp)
        
        xx=xp(nx);%+round(.1*(xp(nx)).^2);
        yy=yp(ny);
        %lattice(yy,xx)=1;
        
        if xx > 3 && xx < dsize(1)-2 && yy > 3 && yy <dsize(2)-2 
            atom=customgauss([ax ax], sx, sy, 0, 0, 1, [dy(yy,xx) dx(yy,xx)]);
            atom=atom/sum(atom(:));
            lattice(yy-2:yy+2,xx-2:xx+2)=atom;
        end
        
    end
end

%lattice=convn(lattice,atom,'same');

end