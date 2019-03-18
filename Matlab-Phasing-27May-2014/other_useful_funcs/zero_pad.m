function [ new_array ] = zero_pad( input,padx,pady,padz )
%does zero padding and cropping
%will pad into 3d if the input is 3d
%padx is for the x direction!  not the first elements
%pady is for the y direction!
%padz is for z, same as usual
%can only do all pad or all crop

try
    padz;
catch
    padz=1;
end

sz=size(input);

type=isreal(input);

dims=ndims(input);

padx=round(padx);
pady=round(pady);
padz=round(padz);

if dims == 2, real_arr=zeros(pady,padx);end
if dims == 3, real_arr=zeros(pady,padx,padz);end

if padx >= sz(2)
    
    padx0=(padx-sz(2));
    pady0=(pady-sz(1));
    
    if dims == 3
        padz0=(padz-sz(3));
        if mod(padz0,2) == 1,ze = 0.5; else ze = 0;end
    end
    
    
    if mod(padx0,2) == 1,xe=0.5;else xe = 0;end         %check for odd numbers
    if mod(pady0,2) == 1,ye=0.5;else ye = 0;end
    
    
%     disp('****************************************')
%     disp('Padding with zeros.....')
%     disp(['original size -' , num2str(sz)])
%     disp(['new size -' , num2str([padx,pady])])
    
    

    real_arr0=real(input); 
    
    if type == 0                    %type =0 means it is complex
        imag_arr0=imag(input);
        
        if dims == 2
            imag_arr=circshift(padarray(imag_arr0,[pady0,padx0],0,'pre'),-[pady0/2+ye,padx0/2+xe]);
            real_arr=circshift(padarray(real_arr0,[pady0,padx0],0,'pre'),-[pady0/2+ye,padx0/2+xe]);
        else
            real_arr=circshift(padarray(real_arr0,[pady0,padx0,padz0],0,'pre'),-[pady0/2+ye,padx0/2+xe,padz0/2+ze]);
            imag_arr=circshift(padarray(image_arr0,[pady0,padx0,padz0],0,'pre'),-[pady0/2+ye,padx0/2+xe,padz0/2+ze]);
        end
        
        new_array=real_arr+1i*imag_arr;
    else
        
        if dims == 2
            real_arr=circshift(padarray(real_arr0,[pady0,padx0],0,'pre'),-[pady0/2+ye,padx0/2+xe]);
        else 
            real_arr=circshift(padarray(real_arr0,[pady0,padx0,padz0],0,'pre'),-[pady0/2+ye,padx0/2+xe,padz0/2+ze]);
        end
        
        new_array=real_arr;
    end
else
    
    %disp('Cropping...')
    
    if mod(sz(2),2) == 1, xc=(sz(2)+1)/2; else xc=sz(2)/2;end  %check for odd
    if mod(sz(1),2) == 1, yc=(sz(1)+1)/2; else yc=sz(1)/2;end  %check for odd
       
    if dims == 3
        if mod(sz(3),2) == 1, zc=(sz(3)+1)/2; else zc=sz(3)/2;end  %check for odd
    end
        
    if mod(padx,2) == 1
        x0=xc-(padx-1)/2;
        x1=(padx-1)/2+xc;
    else
        x0=xc-padx/2+1;
        x1=padx/2+xc;
    end
        
    if mod(pady,2) == 1
        y0=yc-(pady-1)/2;
        y1=(pady-1)/2+yc;
    else
        y0=yc-pady/2+1 ;
        y1=pady/2+yc;
    end    
    
    real_arr0=real(input);
    
    if type == 0                    %check for complex
        
        imag_arr0=imag(input);
        
        if dims == 2
            real_arr=real_arr0(y0:y1,x0:x1);
            imag_arr=imag_arr0(y0:y1,x0:x1);
        else
            z0=sz(3)/2-padz/2+1;
            z1=padz/2+sz(3)/2;
            real_arr=real_arr0(y0:y1,x0:x1,z0:z1);
            imag_arr=imag_arr0(y0:y1,x0:x1,z0:z1);
        end
        new_array = real_arr+1i*imag_arr;
    else
        
        if dims == 2, 
            real_arr=real_arr0(y0:y1,x0:x1); 
        else
                 
            if mod(padz,2) == 1
                z0=zc-(padz-1)/2;
                z1=(padz-1)/2+zc;
            else
                %z0=zc-padz/2;
                %z1=padz/2+zc-1;
                z0=zc-padz/2+1;
                z1=padz/2+zc;
            end  
            real_arr=real_arr0(y0:y1,x0:x1,z0:z1);
            
        end
        
        new_array = real_arr;
        
    end
    
end

real_arr=0;
real_arr0=0;
imag_arr=0;
imag_arr0=0;


    
end

