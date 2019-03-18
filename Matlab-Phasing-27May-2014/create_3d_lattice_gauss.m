function [lattice dx dy dz] = create_3d_lattice_gauss(spacing,dsize,dx,dy,dz,shape)
%jclark
%creates a 2d or 3d latttice
%dsize is x size, ysize,zsize

%gauss 'atoms'
ax=2*spacing+1;
ay=2*spacing+1;
az=2*spacing+1;
sx=.75;
sy=.75;
sz=.75;

hx=floor(ax/2);
hy=floor(ay/2);
hz=floor(az/2);

xp=1:spacing:dsize(1);
yp=1:spacing:dsize(2);
zp=1:spacing:dsize(3);

lattice=zeros([dsize(2),dsize(1),dsize(3)]);

for nx=1:numel(xp)
    for ny=1:numel(yp)
        for nz=1:numel(zp)
            
            xx=xp(nx);
            yy=yp(ny);
            zz=zp(nz);

            if shape(yy,xx,zz) ~= 0
                if xx > floor(ax/2) && xx < dsize(1)-floor(ax/2) && yy > floor(ay/2) && yy <dsize(2)-floor(ay/2) && zz > floor(az/2) && zz <dsize(3)-floor(az/2) 
                    %atom=customgauss([ax ax], sx, sy, 0, 0, 1, [dy(yy,xx) dx(yy,xx)]);
                    atom=gauss_3D_shift(ax,ay,az,sx,sy,sz,dx(yy,xx,zz),dy(yy,xx,zz),dz(yy,xx,zz));
                   
                    %disp([num2str(xx),' ',num2str(yy),' ',num2str(zz)])
                    %disp([num2str(dx(yy,xx,zz)),' ',num2str(dy(yy,xx,zz)),' ',num2str(dz(yy,xx,zz))])
                    
                    %atom=gauss_3D_shift(dsize(1),dsize(2),dsize(3),sx,sy,sz,dx(yy,xx,zz)+xx,dy(yy,xx,zz)+yy,dz(yy,xx,zz)+zz);
                    
                    atom=atom/sum(atom(:));
                    
                    prev=lattice(yy-hy:yy+hy,xx-hx:xx+hx,zz-hz:zz+hz);
                    lattice(yy-hy:yy+hy,xx-hx:xx+hx,zz-hz:zz+hz)=atom+prev;
                    %lattice=lattice+atom;
                    
                    %dx(yy-hy:yy+hy,xx-hx:xx+hx,zz-hz:zz+hz)=dx(yy,xx,zz);
                    %dy(yy-hy:yy+hy,xx-hx:xx+hx,zz-hz:zz+hz)=dy(yy,xx,zz);
                    %dz(yy-hy:yy+hy,xx-hx:xx+hx,zz-hz:zz+hz)=dz(yy,xx,zz);
                    atom=[];
                
                end
                
            end
        end
    end
end

%lattice=convn(lattice,atom,'same');

end