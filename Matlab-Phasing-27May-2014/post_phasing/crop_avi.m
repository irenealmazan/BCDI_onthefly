dir1='/Volumes/HENRIQUE/Users/jesseclark/Documents/MATLAB/data_analysis/DLS-I16-0313-CaCO3/Growth-A/vtks/Pview/';


str1={'1','2','3','4','5'};
str2={'C'};

cd(dir1);

for qq=1:numel(str2)
    for ww=1:numel(str1)
        
        
        string1=['T',char(str1(ww)),'-',char(str2(qq))];
        outstr=[string1,'-Cr.avi'];
        string1=[string1,'.avi'];
        
        cmd=['ffmpeg -i ',string1,' -vf crop=400:400:30:50 -vcodec mjpeg -qscale 2 ',outstr];
        
        system(cmd)
        
    end
    
end



1;
