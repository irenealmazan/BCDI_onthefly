function [ params ] = determine_resolution( data,dnc )
%jclark
%dnc is the amount from center


params.dsx=1;params.dsy=1;params.dsz=1;

orient={'xy','xz','zy'};

if exist('dnc') ~= 1
    params.xy=extract_3D_slice(data,'xy');
    params.xz=extract_3D_slice(data,'xz');
    params.zy=extract_3D_slice(data,'zy');
else  %do them off center if dnc is specified
    params.xy=extract_3D_slice(data,'xy',ceil(size(data,3)/2)+dnc);
    params.xz=extract_3D_slice(data,'xz',ceil(size(data,1)/2)+dnc);
    params.zy=extract_3D_slice(data,'zy',ceil(size(data,2)/2)+dnc);
end


line_avg=5;%ensure this is odd

params.xline=mean(params.xy(size(params.xy,1)/2-(line_avg-1)/2:size(params.xy,1)/2+(line_avg-1)/2,:),1);

params.yline=mean(params.xy(:,size(params.xy,2)/2-(line_avg-1)/2:size(params.xy,2)/2+(line_avg-1)/2),2);

params.zline=mean(params.zy(size(params.zy,1)/2-(line_avg-1)/2:size(params.zy,1)/2+(line_avg-1)/2,:),1);

params.dx=diff(params.xline);

params.dy=diff(params.yline);

params.dz=diff(params.zline);

%extract just the edges
width=8;    %number of pixels to fit the gaussian to
params.dx_e1=params.dx( find(params.dx == max(params.dx))-width/2:find(params.dx == max(params.dx))+width/2);
params.dx_e2=-params.dx( find(params.dx == min(params.dx))-width/2:find(params.dx == min(params.dx))+width/2);

params.dy_e1=params.dy( find(params.dy == max(params.dy))-width/2:find(params.dy == max(params.dy))+width/2);
params.dy_e2=-params.dy( find(params.dy == min(params.dy))-width/2:find(params.dy == min(params.dy))+width/2);

params.dz_e1=params.dz( find(params.dz == max(params.dz))-width/2:find(params.dz == max(params.dz))+width/2);
params.dz_e2=-params.dz( find(params.dz == min(params.dz))-width/2:find(params.dz == min(params.dz))+width/2);

[params.fit_dx_e1 params.parms_xe1]=fit_gauss_data(params.dx_e1);
[params.fit_dx_e2 params.parms_xe2]=fit_gauss_data(params.dx_e2);

[params.fit_dy_e1 params.parms_ye1]=fit_gauss_data(params.dy_e1);
[params.fit_dy_e2 params.parms_ye2]=fit_gauss_data(params.dy_e2);

[params.fit_dz_e1 params.parms_ze1]=fit_gauss_data(params.dz_e1);
[params.fit_dz_e2 params.parms_ze2]=fit_gauss_data(params.dz_e2);



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