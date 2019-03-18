function [ params ] = determine_resolution_2D(xline,yline)
%jclark

params.dsx=1;params.dsy=1;params.dsz=1;



params.xline=xline;

params.yline=yline;

params.dx=diff(params.xline);

params.dy=diff(params.yline);



%extract just the edges
width=8;    %number of pixels to fit the gaussian to

ll=find(params.dx == max(params.dx))-width/2;
rr=find(params.dx == max(params.dx))+width/2;
if ll < 1,ll=1;end
if rr > max(size(params.dx)),rr=max(size(params.dx));end
params.dx_e1=params.dx( ll:rr);

ll=find(params.dx == min(params.dx))-width/2;
rr=find(params.dx == min(params.dx))+width/2;
if ll < 1,ll=1;end
if rr > max(size(params.dx)),rr=max(size(params.dx));end
params.dx_e2=-params.dx( ll:rr);

ll=find(params.dy == max(params.dy))-width/2;
rr=find(params.dy == max(params.dy))+width/2;
if ll < 1,ll=1;end
if rr > max(size(params.dy)),rr=max(size(params.dy));end
params.dy_e1=params.dy( ll:rr);

ll=find(params.dy == min(params.dy))-width/2;
rr=find(params.dy == min(params.dy))+width/2;
if ll < 1,ll=1;end
if rr > max(size(params.dy)),rr=max(size(params.dy));end
params.dy_e2=-params.dy( ll:rr);


[params.fit_dx_e1 params.parms_xe1]=fit_gauss_data(params.dx_e1);
[params.fit_dx_e2 params.parms_xe2]=fit_gauss_data(params.dx_e2);

[params.fit_dy_e1 params.parms_ye1]=fit_gauss_data(params.dy_e1);
[params.fit_dy_e2 params.parms_ye2]=fit_gauss_data(params.dy_e2);

params.sig_x=mean([params.parms_xe1(2),params.parms_xe2(2)]);
params.sig_y=mean([params.parms_ye1(2),params.parms_ye2(2)]);


end
function [gg x] = fit_gauss_data(data)

options = optimset('Display','off','Algorithm','interior-point');
 
x0(1)=find(data == max(data));

x0(2)=2;

x0(3)=max(data);

%x0(4)=-1/5;
%x0(5)=1/20*x0(3);

lb=[0,0,0];%,-5,-x0(3)];
ub=[numel(data),numel(data),5*x0(3)];%,5,x0(3)];

f=@(x)gauss_fit(x,data);
x=fmincon(f,x0,[],[],[],[],[lb],[ub],[],options);


mean=x(1);
sig=x(2);
A0=x(3);
%m=x(4);
%c=x(5);

xx=1:numel(data);   %abisca values


gg=A0*exp(-0.5*(xx-mean).^2/sig.^2);%+m*xx+c;


end
function E = gauss_fit(x,data,range )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%y=A0*exp(-0.5 x*x/(sig*sig))+mx+c

mean=x(1);
sig=x(2);
A0=x(3);
%m=x(4);
%c=x(5);


xx=1:numel(data);   %abisca values


gauss=A0*exp(-0.5*(xx-mean).^2/sig.^2);%+m*xx+c;

data=data(:);
gauss=gauss(:);

E=sum(abs(gauss-data));


end