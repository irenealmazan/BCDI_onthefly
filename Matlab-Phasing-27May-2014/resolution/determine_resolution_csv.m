function determine_resolution_csv()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


dir1='/Users/jesseclark/Documents/MATLAB/data_analysis/Au0710/Crystal1/261-50-2/rand-starts/Rec-Rnd3-261-ERHIO200-CVl-SW/';

colors={'blue','red'};

dir_file=dir1;
pfile=rdir([dir_file,'**/*.csv']);
sz=size(pfile);

lw=1.5;
font_size=25;
clr=char(colors(1));

norm=0;

save_dir=dir1;

for qq=1:sz(1)

        dname0=pfile(qq).name;
        dname=dname0(end-7:end);

        disp(dname)

        switch dname

            case 'x_lo.csv'
                xd=csvread(dname0,1,0);
                x_scalars=xd(:,1);
                x_xx=xd(:,4);
                
                ind=logical(1-isnan(x_scalars));
                
                x_s=x_scalars(ind);
                x_x=x_xx(ind);
                x_scalars=x_s;
                x_xx=x_x;
                
                if norm == 1,x_scalars=x_scalars/max(x_scalars);end
                
                qx=x_xx;
                x=x_scalars;
                dx=abs(qx(2)-qx(1));

            case 'y_lo.csv'
                yd=csvread(dname0,1,0);
                y_scalars=yd(:,1);
                y_yy=yd(:,4);
                
                ind=logical(1-isnan(y_scalars));
                
                y_s=y_scalars(ind);
                y_y=y_yy(ind);
                y_scalars=y_s;
                y_yy=y_y;
                
                if norm == 1,y_scalars=y_scalars/max(y_scalars);end
                
                qy=y_yy;
                y=y_scalars;
                dy=abs(qy(2)-qy(1));
              

            case 'z_lo.csv'
                zd=csvread(dname0,1,0);
                z_scalars=zd(:,1);
                z_zz=zd(:,4);
                
                ind=logical(1-isnan(z_scalars));
                
                z_s=z_scalars(ind);
                z_z=z_zz(ind);
                z_scalars=z_s;
                z_zz=z_z;
                
                if norm == 1,z_scalars=z_scalars/max(z_scalars);end
                
                qz=z_zz;
                z=z_scalars;
                dz=abs(qz(2)-qz(1));
                
        end
                
end


params.xline=x_scalars;
params.yline=y_scalars;
params.zline=z_scalars;

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
                
output_resolution(params,[save_dir,'csv-'])                

end

