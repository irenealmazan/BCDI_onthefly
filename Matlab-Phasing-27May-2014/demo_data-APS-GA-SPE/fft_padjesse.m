function data = fft_padjesse(data,bin)

%% set an array of values that it will automtically pad to to give 
%a good number for the FFT, i.e powers of 2 or PX2^N, where P is a low
%numbered prime i.e <7 or 9
%pref=[32,48,64,80,96,128,144,160,192,256,320,512];
pow_2=2.^(1:10);
pref=sort([2.^(4:16),2.^(3:10)*3,2.^(3:10)*5,2.^(3:10)*9]);%sort([2.^(6:10),pow_2*3,pow_2*5,pow_2*9]);

dsz=size(data);

disp('')

try
    bin;
catch
    bin=[1,1,1];
end

if numel(dsz) == 3
    sx=dsz(2);
    sy=dsz(1);
    sz=dsz(3);

    disp(['Before padding for FFT [x,y,z] - [',num2str([sx,sy,sz]),']'])
    %disp(['Array size after binning [x,y,z] - [',num2str([sx,sy,sz]),']'] )

    while sum(sx == bin(1)*pref) == 0,sx=sx+1;end
    while sum(sy == bin(2)*pref) == 0,sy=sy+1;end
    while sum(sz == bin(3)*pref) == 0,sz=sz+1;end

    nn=[sx,sy,sz];
    disp(['Padding for FFT [x,y,z] - [',num2str([sx,sy,sz]),']'])
    %disp(['After binning [x,y,z] - [',num2str([sx/bin(1),sy/bin(2),sz/bin(3)]),']'])

    %data=zero_pad_ver2(data,nn(1),nn(2),nn(3));
    data=zero_pad_ver3(data,nn(1),nn(2),nn(3));
end
if numel(dsz) == 2
    sx=dsz(2);
    sy=dsz(1);

    disp(['Before padding for FFT [x,y] - [',num2str([sx,sy]),']'])
    %disp(['Array size after binning [x,y] - [',num2str([sx,sy]),']'] )

    while sum(sx == bin(1)*pref) == 0,sx=sx+1;end
    while sum(sy == bin(2)*pref) == 0,sy=sy+1;end

    nn=[sx,sy];
    disp(['Padding for FFT [x,y] - [',num2str([sx,sy]),']'])
    %disp(['After binning [x,y] - [',num2str([sx/bin(1),sy/bin(2)]),']'])


    data=zero_pad_ver2(data,nn(1),nn(2));
end

disp('Complete.....')
disp(' ')

%disp('----------------------------------------------------------------')


end
