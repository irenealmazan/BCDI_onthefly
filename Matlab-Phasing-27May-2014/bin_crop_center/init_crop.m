function data = init_crop(data,nnc) 
%jclark

%% initial cropping

if numel(nnc) == 6
    
    
    %if sum(nnc(1:2)) > 0, data=crop_dim(data,nnc(1:2),1);end
    %if sum(nnc(3:4)) > 0, data=crop_dim(data,nnc(3:4),2);end
    %if sum(nnc(5:6)) > 0, data=crop_dim(data,nnc(5:6),3);end
    
    sz=size(data);
    xs=sz(2);
    ys=sz(1);

    if ndims(data) == 3,zs=sz(3);end
    
    if sum(nnc(1:2)) < 0,
        %disp('doing intial x cropping....')
        nncc=0;
        nncc=[abs(nnc(1))+1,xs-abs(nnc(2))];
       % disp(['xs [',num2str(xs),'] --> [',num2str(1+nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,1);end
    
    if sum(nnc(3:4)) < 0, 
        %disp('doing intial y cropping....')
        nncc=0;
        nncc=[abs(nnc(3))+1,ys-abs(nnc(4))];
     %   disp(['ys [',num2str(ys),'] --> [',num2str(1+nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,2);end
   
    
    if sum(nnc(5:6)) < 0, 
        %disp('doing intial z cropping....')
        nncc=0;
        nncc=[abs(nnc(5))+1,zs-abs(nnc(6))];
      %  disp(['zs [',num2str(zs),'] --> [',num2str(1+nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,3);end
else
    disp('if initial cropping is required')
    disp('set nnc=[-x0,-x1,-y0,-y1,-z0,-z1], where nnc gives pixels to crop')
    disp('off each dimension otherwise set nnc=[0]')
end

end