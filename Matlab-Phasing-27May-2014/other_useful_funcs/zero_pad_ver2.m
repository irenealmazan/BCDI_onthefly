function [ new_array ] = zero_pad_ver2( input,padx,pady,padz )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

sz=size(input);

type=isreal(input);

dims=ndims(input);


real_arr=isreal(input);

i=complex(0,1);

if dims == 2,real_array=all_pad(real(input),padx,pady);end
if dims == 3,real_array=all_pad(real(input),padx,pady,padz);end

if real_arr== 0
    if dims == 2,imag_array=all_pad(imag(input),padx,pady);end
    if dims == 3,imag_array=all_pad(imag(input),padx,pady,padz);end
    new_array=real_array+i*imag_array;
else new_array=real_array;end
    

    
function [new_arr] = all_pad(input,padx,pady,padz)
    
    sn=size(input);
    dim=ndims(input);
    
    if dim == 3,padn=[padx,pady,padz]-[sn(2),sn(1),sn(3)];end  %the amount to pad 
    if dim == 2,padn=[padx,pady]-[sn(2),sn(1)];end
    
    pre_pad=zeros(1,dim);
    pos_pad=zeros(1,dim);
    
    for qq = 1:dim
        
        odd=mod(padn(qq),2);
        
        if odd == 1, 
            pre_pad(qq)=(padn(qq)-1)/2;
            pos_pad(qq)=(padn(qq)+1)/2;
        else
            pre_pad(qq)=padn(qq)/2;
            pos_pad(qq)=padn(qq)/2;
        end
        
    end
    
    temp_pad=pre_pad;
    pre_pad(1)=temp_pad(2);
    pre_pad(2)=temp_pad(1);
    temp_pad=pos_pad;
    pos_pad(1)=temp_pad(2);
    pos_pad(2)=temp_pad(1);
    
    new_arr=padarray(input,pre_pad,0,'pre');
    new_arr=padarray(new_arr,pos_pad,0,'post');
        
end


end

