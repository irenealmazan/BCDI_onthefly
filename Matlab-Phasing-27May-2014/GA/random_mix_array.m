function [array3 array4] = random_mix_array( array1,array2,crossover )
%jclark
%randomly mixes elements of arrays

try
    crossover;
catch
    crossover=0.5;
end


ind1=array1(:); %get arr1 vals
ind2=array2(:); %same for 2

cind=random('uniform',0,1,size(ind1));   %get random indices

cind( cind > crossover)=0;
cind( cind ~= 0)=1;

cind=logical(cind);

array3=array1;
array4=array2;

array3(cind)=array2(cind);
array4(cind)=array1(cind);

array1=[];
array2=[];

end

