function data = init_pad(data,nnc)
%jclark

sz=size(data);
xs=sz(2);
ys=sz(1);

try
 zs=sz(3);
end
 
if numel(nnc) == 6
    
    if sum(nnc(1:2)) > 0, 
        %disp('doing intial x padding....')
        data=padarray(data,[0 abs(nnc(1)) 0],0,'pre');
        data=padarray(data,[0 abs(nnc(2)) 0],0,'post');
        %disp(['xs [',num2str(xs),'] --> [',num2str(xs+nnc(2)+nnc(1)),']'] )
    end
    
    if sum(nnc(3:4)) > 0, 
        %disp('doing intial y padding....')
        data=padarray(data,[abs(nnc(3)) 0 0],0,'pre');
        data=padarray(data,[abs(nnc(4)) 0 0],0,'post');
        %disp(['ys [',num2str(ys),'] --> [',num2str(ys+nnc(4)+nnc(3)),']'] )
    end
    
    if sum(nnc(5:6)) > 0, 
        %disp('doing intial z padding....')
        data=padarray(data,[0 0 abs(nnc(5))],0,'pre');
        data=padarray(data,[0 0 abs(nnc(6))],0,'post');
        %disp(['zs [',num2str(zs),'] --> [',num2str(zs+nnc(5)+nnc(6)),']'] )
    end
    
else
    disp('if initial pading is required')
    disp('set nnc=[+x0,+x1,+y0,+y1,+z0,+z1], where nnc gives pixels to pad')
    disp('onto each dimension otherwise set nnc=[0]')
end


end
