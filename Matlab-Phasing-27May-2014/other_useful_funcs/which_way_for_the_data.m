function [ out_flip ] = which_way_for_the_data(data1,data2)
%jclark
%data1 is reference
%data2 is to be aligned

qq=0;


for dim1=1:2
    for dim2=1:2
        for dim3=1:2
            
            dtemp=data2;
            qq=qq+1;
            
            if dim1==1,dtemp=flipdim(dtemp,1);flip1=1;else flip1=0;end
            if dim2==1,dtemp=flipdim(dtemp,2);flip2=1;else flip2=0;end
            if dim3==1,dtemp=flipdim(dtemp,3);flip3=1;else flip3=0;end
            
            flip(qq,:)=[flip1,flip2,flip3];
            
            nxcc(qq)=max(max(max(normxcorr3(dtemp, data1, 'same'))));
            dtemp=[];
        end
    end
end

ind=(nxcc == max(nxcc(:)));
out_flip=flip(ind,:);

end

