function [ temp ] = ResizeFFt( first,temp )
%give two arrays, will resize in the fft domain rather than realspace
%first is the reference and temp is the one to resize

if sum(size(temp)-size(first)) ~= 0
    
    %disp('Resizing for comparison....')
    tempFFt=fftshift(fftn(fftshift(temp)));
    tempFFt=zero_pad_ver3(tempFFt,size(first,2),size(first,1),size(first,3));
    temp=fftshift(ifftn(fftshift(tempFFt)));
                
end


end

