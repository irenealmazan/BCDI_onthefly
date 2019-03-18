dir1='/Volumes/HENRIQUE/Users/jesseclark/Documents/MATLAB/data_analysis/DLS-I16-0313-CaCO3/Growth-A/vtks/Pview/';

%dir11=make_string_from_recursions([dir1,'T'],{'1-A','1-B','1-C'},'/A*');

%dir22=make_string_from_recursions([dir1,'T'],{'2-A','2-B','2-C'},'/A*');

%dir33=make_string_from_recursions([dir1,'T'],{'3-A','3-B','3-C'},'/A*');

str1={'1','2','3','4','5'};
str2={'C'};

for qq=1:numel(str2)
    for ww=1:numel(str1)
        
        
        string1=['T',char(str1(ww)),'-',char(str2(qq))];
        string2=[dir1,string1,'/A*'];
        
        make_avi_from_images(string2,[dir1,string1,'.avi'],30,500);
    end
    
end



1;
