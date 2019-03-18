function [ self ] = UpdateCoordSystem(self,T)
%%jclark
%self.size1
%self.size2
%self.size3
%
%%
%def UpdateCoordSystem(self, dims):

%print "dims in getcoord", dims
%if len(dims) < 2:
%    return

%self.size1 = dims[0]

%self.size2 = dims[1]

%self.size3 = dims[2]
%%
Nx=self.size1;
Ny=self.size2;
Nz=self.size3;
self.dx=1./self.size1;
self.dy=1./self.size2;
self.dz=1./self.size3;

%r=nd.mgrid[ (self.size1-1)*self.dx:-self.dx:-self.dx, ...
%        0:self.size2*self.dy:self.dy, 0:self.size3*self.dz:self.dz];

[x y z]=meshgrid( ((1:self.size1)*self.dx),((1:self.size2)*self.dy),(1:self.size3)*self.dz);
%y=flipdim(y,2);
%x=flipdim(x,1);
%[x y z]=meshgrid( reverse((1:self.size1)*self.dx),reverse((1:self.size2)*self.dy),(1:self.size3)*self.dz);
%[x y z]=meshgrid( (self.size1-1):-self.dx,(1:self.size2)*self.dy,(1:self.size3)*self.dz);

self.x=x;
self.y=y;
self.z=z;

r=zeros([3,Ny,Nx,Nz]);

r(1,:,:,:)=x;
r(2,:,:,:)=y;
r(3,:,:,:)=z;

%r=zeros([3,Nx,Ny,Nz]);
%r(1,:,:,:)=x;
%r(2,:,:,:)=y;
%r(3,:,:,:)=z;
%r=reshape(r,3,self.size1*self.size2*self.size3);

r=reshape(r,3,self.size2*self.size1*self.size3);
%r.shape=3,self.size1*self.size2*self.size3

%r=r.transpose()
%r = transpose(r);

%self.coords = nd.dot(r, self.T)
self.coords=T*r;

self.coords=reshape(self.coords,3,self.size2,self.size1,self.size3);
%return self.coords



end

