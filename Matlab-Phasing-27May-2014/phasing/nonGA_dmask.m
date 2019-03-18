function gmask = nonGA_dmask(params,sigxy)
%jclark

if sigxy < 1,
    disp(['USING GAUSS MASK - ',num2str(sigxy)])
    if params.do_2D == 1,
        gx=params.nn(1);
        gy=params.nn(2);
        gmask=gauss_2D(gx,gy,gx*sigxy,gy*sigxy);
        gmask=gmask/max(gmask(:));
        %gmask=sinc_2D(gx,gy,sigxy*2,sigxy*2);
        %gmask=super_gauss_2D(gx,gy,gx*sigxy,gy*sigxy,2);
    else
        gx=params.nn(1);
        gy=params.nn(2);
        gz=params.nn(3);
        gmask=gauss_3D(gx,gy,gz,gx*sigxy,gy*sigxy,gz*sigxy);
        gmask=gmask/max(gmask(:));
    end
else
    gmask=1; 
end


end