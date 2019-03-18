function [support_i] = zero_pad_ver_Irene(arrysize,arrysize2,arrysize3,data_to_pad)
% this function turns an array/matrix of a certain size, to a zero padded
% and centered array/matrix of size size_to_pad

     support_1 = data_to_pad;
     support_i = zeros(arrysize,arrysize2,arrysize3);
     
    
     support_i(round(arrysize/2-size(support_1,1)/2)+1:round(arrysize/2+size(support_1,1)/2),...
         round(arrysize2/2-size(support_1,1)/2)+1:round(arrysize2/2+size(support_1,1)/2), ...
         round(arrysize3/2-size(support_1,3)/2)+1:round(arrysize3/2+size(support_1,3)/2))...
         = support_1;


   
    
end