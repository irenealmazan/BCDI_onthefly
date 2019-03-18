function [Fxy] = extract_defect_netwrok(pnm)
%jclark

nangles=0:1:360;

pnm=pnm.*shrink_wrap(abs(pnm),.1,1);

for ww=1:numel(nangles)

    ph00=exp(i*deg2rad(nangles(ww)));
    
    %calc gradient
    [Ax,Ay,Az,Fx,Fy,Fz,ss] = calc_gradient(abs(pnm),angle(ph00.*pnm),2);
    Ax=Ax.*ss;Ay=Ay.*ss;Az=Az.*ss;
    Fx=Fx.*ss;Fy=Fy.*ss;Fz=Fz.*ss;

    %calc gradient on gradient
    [Axx,Ayx,Azx] = calc_gradient((Fx));
    [Axy,Ayy,Azy] = calc_gradient((Fy));
    [Axz,Ayz,Azz] = calc_gradient((Fz));

    A=abs(Axx.*Ayx.*Azx);
    B=abs(Axy.*Ayy.*Azy);
    C=abs(Axz.*Ayz.*Azz);
    D=max(cat(4,A,B,C),[],4);

    %Threshold
    D=D./max(D(:));
    Ds=std(D(D ~= 0));
    D(D < Ds)=0;

    D=abs(convn(D,gauss_3D(11,11,11,.75,.75,.75),'same'));
    E=abs(convn(A+B+C,gauss_3D(11,11,11,.75,.75,.75),'same'));

    Fxy=D;%D
    Fxy=Fxy./max(abs(Fxy(:)));

    %get the connected regions
    BW=(Fxy >= .02);
    CC = bwconncomp(BW,18);

    %keep ones that are above 10 pixels
    for ll=1:numel(CC.PixelIdxList)
        if numel(CC.PixelIdxList{ll}) <  25,BW(CC.PixelIdxList{ll})=0;end
    end
    Fxy=Fxy.*BW;

    if ww == 1,Nxy=Fxy;else
        Nxy=cat(4,Nxy,Fxy);
    end
end
    
    
    Fxy=min(Nxy,[],4);
    
    Fxy=Fxy./(max(Fxy(:)));
    
    BW=(Fxy >= .01);
    CC = bwconncomp(BW,18);
    for ll=1:numel(CC.PixelIdxList)
            if numel(CC.PixelIdxList{ll}) <  25,BW(CC.PixelIdxList{ll})=0;end
        end
    Fxy=Fxy.*BW;

end

