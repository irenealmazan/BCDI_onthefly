function [ up_arr ] = FFT_interpolate(arr,new_n )
%% interpolates using an FFT
% arr is the in array
% new_n is the new size in pixels (should be even)
% new_n=[new_x,new_y,new_z];  leave new_z out to do 2d

%%
nn=size(arr);
sz=nn;

sz(1)=nn(2);
sz(2)=nn(1);

if ndims(arr) == 2
    
    pad=new_n-sz;   %determine amount to pad
    if mod(pad(1),2) == 1; pad(1)=pad(1)+1;end   %make it even if it is odd
    if mod(pad(2),2) == 1; pad(2)=pad(2)+1;end   %make it even if it is odd
    
    new_n=pad+sz;
    pad=pad/2;
    
    nnc=[pad(1),pad(1),pad(2),pad(2)];
    
    fft_arr=fftshift(fftn(ifftshift(arr)));
    
    fft_arr=init_pad(fft_arr,nnc);
    fft_arr=init_crop(fft_arr,nnc);
    
    up_arr=fftshift(ifftn(ifftshift(fft_arr)));
    
end

if ndims(arr) == 3
    
    pad=new_n-sz;   %determine amount to pad
    if mod(pad(1),2) == 1; pad(1)=pad(1)+1;end   %make it even if it is odd
    if mod(pad(2),2) == 1; pad(2)=pad(2)+1;end   %make it even if it is odd
    if mod(pad(3),2) == 1; pad(3)=pad(3)+1;end   %make it even if it is odd
    
    new_n=pad+sz;
    pad=pad/2;
    
    nnc=[pad(1),pad(1),pad(2),pad(2),pad(3),pad(3)];
    
    fft_arr=fftshift(fftn(ifftshift(arr)));
    
    fft_arr=init_pad(fft_arr,nnc);
    fft_arr=init_crop(fft_arr,nnc);
    
    up_arr=fftshift(ifftn(ifftshift(fft_arr)));
    
end

upscale=new_n/sz;

disp(['Upscaling factor - [',num2str(upscale),']'])

end

function data = init_pad(data,nnc)

if ndims(data) == 3
    sz=size(data);
    xs=sz(2);
    ys=sz(1);
    zs=sz(3);

    if numel(nnc) == 6

        if sum(nnc(1:2)) > 0, 
            disp('doing intial x padding....')
            data=padarray(data,[0 abs(nnc(1)) 0],0,'pre');
            data=padarray(data,[0 abs(nnc(2)) 0],0,'post');
            disp(['xs [',num2str(xs),'] --> [',num2str(xs+nnc(2)+nnc(1)),']'] )
        end

        if sum(nnc(3:4)) > 0, 
            disp('doing intial y padding....')
            data=padarray(data,[abs(nnc(3)) 0 0],0,'pre');
            data=padarray(data,[abs(nnc(4)) 0 0],0,'post');
            disp(['ys [',num2str(ys),'] --> [',num2str(xs+nnc(4)+nnc(3)),']'] )
        end

        if sum(nnc(5:6)) > 0, 
            disp('doing intial z padding....')
            data=padarray(data,[0 0 abs(nnc(5))],0,'pre');
            data=padarray(data,[0 0 abs(nnc(6))],0,'post');
            disp(['zs [',num2str(zs),'] --> [',num2str(zs+nnc(5)+nnc(6)),']'] )
        end

    end
end

if ndims(data) == 2
    sz=size(data);
    xs=sz(2);
    ys=sz(1);
    

    if numel(nnc) >= 4

        if sum(nnc(1:2)) > 0, 
            disp('doing intial x padding....')
            data=padarray(data,[0 abs(nnc(1)) 0],0,'pre');
            data=padarray(data,[0 abs(nnc(2)) 0],0,'post');
            disp(['xs [',num2str(xs),'] --> [',num2str(xs+nnc(2)+nnc(1)),']'] )
        end

        if sum(nnc(3:4)) > 0, 
            disp('doing intial y padding....')
            data=padarray(data,[abs(nnc(3)) 0 0],0,'pre');
            data=padarray(data,[abs(nnc(4)) 0 0],0,'post');
            disp(['ys [',num2str(ys),'] --> [',num2str(xs+nnc(4)+nnc(3)),']'] )
        end

    end
end

end