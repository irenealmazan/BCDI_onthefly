function [ output_args ] = hsv_color_whel( input_args )
%jclark


nn=512;
rad=400;

[x y]=meshgrid(-nn:nn,-nn:nn);

r=sqrt(x.^2+y.^2);

ind=(r < rad);

tr=r.*ind;

tr(r > rad)=0;

tr=tr-min(tr(:));

tr(r > rad)=0;

phase=atan2(y,x);

phase(r >rad)=0;

tr=tr.*exp(-i*phase);



end

