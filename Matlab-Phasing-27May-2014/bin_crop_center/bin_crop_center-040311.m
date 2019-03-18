function [ data ] = bin_crop_center(files,bgs,bin,mindata,aliens,nnc,file_params)
% - Jesse Clark, LCN, UCL October-November 2010
%   jesse.clark@ucl.ac.uk, jessenclark@gmail.com
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%Can read SPE, tif and mat files.
%
%For .spe files it is assumed that it is a 3d data set.  enter file names
%into files files={blah.spe,blablah.spe}.  background files are the same.
%doesn't matter if these are 2d or 3d, will just choose one background
%image for subtraction.  leave empty for no bg subtraction.
%
%for .tif files, it is assumed these have come from a pilatus and the 3d
%data set is made up of n 2d images.  no bg 
%subtraction is done.  just enter the first and last file names and the
%number of accuulations per theat position (file_params.nfiles=#) and it
%will make the 3d data set.
%
%alternatively provide a .mat file and it will directly read in a matlab
%file, assumed to be already in 3d.
%
%providing .spe,.tif or.mat will still use mindata do remove values below
%this and will center the data onto the maximum and remove aliens.

disp('----------------------------------------------------------------')
%disp('Loading files....')

if isempty(strfind(lower(char(files{1})),'.spe')) == 0
    disp('Loading .SPE data files....')
    
    nfiles=size(files);             %get the number of files
    nfiles=max(nfiles);

    disp([num2str(nfiles),' will be loaded...'])

    disp(['Values below -[',num2str(mindata),']- will be set to 0...'])
    disp(' ')
    disp('Loading file - [1]')
    data0=speread(char(files(1)));  %read first file in

    data0=data0.data;               %keep data only

    nz=size(data0);

    if numel(nz) ~= 2,nz=nz(3);end

    data0=rot3d(data0);

    nx=size(data0);                 %get dimensions

    if numel(bgs) == 0,back0=0;else   %check if a bg files is supplied
        back0=speread(char(bgs(1)));          %get background image
        disp('Loading background image - [1] ')
        back0=rot3d(back0.data);

        if ndims(back0) == 3,back0=back0(:,:,2);end
        %if exist('back_no') == 0,back_no=ndims(back0)-1;end
        %if back_no == 2,back0=back0(:,:,back_no);end
        %if back_no == 1,back0=back0(:,:,back_no);end
    end     
    %disp('Subtracting background image from data...')

    if numel(nx) == 2,nx=[nx,1];end

    data=data0-data0;
    for qq = 1:nx(3), data(:,:,qq)=flipud(data0(:,:,qq)-back0);end         %remove bg
    data0=0;
    back0=0;

    ind=(data < mindata);             %keep everything above mindata
    data(ind)=0;

    %get center of 1st data set so others can be shifted to a common center
    %for alien removal
    [u1 u2 u3]=ind2sub(size(data),find(data == max(max(max(data)))));  %get center of 1st data set
    if numel(u1) ~= 1,
        saturation_warning(char(files(ww)));
        u1=u1(1);
        u2=u2(1);
        u3=u3(1);
    end
    
    if nfiles > 1                   % load other files if present

        for ww = 2:nfiles
            disp(['Loading file - [',num2str(ww),']'])
            data1=speread(char(files(ww)));
            data1=rot3d(data1.data);

            if numel(bgs) == 0,back1=0;else
                disp(['Loading background image - [',num2str(ww),']'])
                back1=speread(char(bgs(ww)));          %get background image
                back1=rot3d(back1.data);

                %if exist('back_no') == 0,back_no=ndims(back1)-1;end
                %if back_no == 2,back1=back1(:,:,back_no);end

                if ndims(back1) == 3,back1=back1(:,:,2);end
                %back1=back1(:,:,2);
            end
            
            for qq = 1:nx(3), data1(:,:,qq)=flipud(data1(:,:,qq)-back1);end  
           
            ind=( data1 < mindata );
            data1(ind)=0;
            %data1=center_array(data1);
            
            [v1 v2 v3]=ind2sub(size(data1),find(data1 == max(max(max(data1)))));  %get center of next data set
            
            if numel(v1) ~= 1,
                saturation_warning(char(files(ww)));
                v1=v1(1);
                v2=v2(1);
                v3=v3(1);
            end
            
            data1=circshift(data1,[u1-v1,u2-v2,u3-v3]);
            
            data=data+data1;
        end

    end
    
    data=remove_aliens(data,aliens);
    data = init_pad(data,nnc);
    data = init_crop(data,nnc);
    data=fft_pad(data,[bin,1]);
    data=center_array(data);
end

if isempty(strfind(lower(char(files{1})),'.tif')) == 0
    disp('Loading .tif data files....')
    idx = regexp(char(files{1}),'\d+');
    nums = regexp(char(files{1}),'\d+','match');   %get the numbers from the file name
    
    first=char(nums(end));        %get first 
    
    idx = regexp(char(files{2}),'\d+');
    nums = regexp(char(files{2}),'\d+','match');   %get the numbers from the file name
    
    last=char(nums(end));
    
    file0=char(files{1});                   %get the prefix to generate other file names
    prefix=file0(1:length(file0)-length(char(nums(end)))-length('.tif'));
    
    last=str2num(char(last));       %convert it to a number
    first=str2num(char(first));
    total_files=last-first+1;       %number of files that will be loaded
    
    disp([num2str(total_files),' will be loaded...'])

    disp(['Values below -[',num2str(mindata),']- will be set to 0...'])
    disp(' ')
    
    temp = imread(file0, 'tif');
    
    ntheta=total_files/file_params.nfiles;
    
    data=zeros([size(temp),ntheta]);    %create data array
    nx=size(data);
    temp=0;
    counter=0;
    
    for qq = 1:ntheta
        
        temp0=0;
        
        for pp = 1:file_params.nfiles
 
            fname=[prefix,num2str(first+counter),'.tif'];
            disp(fname(end-10:end))
            
            temp = imread(fname, 'tif');
            
            temp0=temp0+temp;
            temp=0;
            counter=counter+1;
        end
        
        data(:,:,qq)=temp0;
    
    end
    
    data=remove_aliens(data,aliens);
    data = init_pad(data,nnc);
    data = init_crop(data,nnc);
    data=fft_pad(data,[bin,1]);
    data=center_array(data);
    
    
    ind=( data < mindata );
    data(ind)=0;
end

if isempty(strfind(lower(char(files{1})),'.mat')) == 0
    
    disp('Loading MATLAB data file....')
    
    load(char(files{1}));
    %data=center_array(data);
    
    nx=size(data);
    disp(['Values below -[',num2str(mindata),']- will be set to 0...'])
    disp(' ')
    ind=( data < mindata );
    data(ind)=0;
    
    data=remove_aliens(data,aliens);
    data = init_pad(data,nnc);
    data = init_crop(data,nnc);
    data=fft_pad(data,[bin,1]);
    data=center_array(data);
    
    
end


%%

nx=size(data);
%%
disp(' ')
figure(89)
disp('Plotting histogram before binning....')
ind=(data <2000);
[h x]=hist(data(ind),100);
%p1=plot(x,h);%
p1=plot(x,log(h));
ymin=log(h(2:end));
ymin=min(ymin(ymin>0));
axis([0 2000 ymin max(log(h(2:end)))])
set(p1,'Color','blue','LineWidth',1.0);
xlabel('ADU''s');
ylabel('frequency');
title('data histogram');

%pad the array so that it bins exactly

x0=nx(2);
y0=nx(1);

disp(' ')
disp('Resizing data....')
disp(['Original data size [x,y,z] - [',num2str([x0,y0,nx(3)]),']'])

while mod(x0,bin(1)) ~=0,x0=x0+1;end
while mod(y0,bin(2)) ~=0,y0=y0+1;end


data=padarray(data,[(y0-nx(1)),(x0-nx(2)),0],0,'pre');

%data=imresize(center_array(data),[y0/bin(2),x0/bin(1)],'nearest');
data_new=zeros([y0/bin(2),x0/bin(1),nx(3)]);
for qq = 1:nx(3), 
    
    data_new(:,:,qq)=box_interp(data(:,:,qq),bin(1),bin(2));
end
data=data_new;
data_new=0;

disp(['Array size after binning [x,y,z] - [',num2str([x0/bin(1),y0/bin(2),nx(3)]),']'])

try
    disp(['doing secondary thresholding - [',num2str(file_params.schot_th),']'])
    ind=(data < file_params.schot_th);
    data(ind)=0;
end

disp(' ')
figure(90)
disp('Plotting histogram after binning....')
ind=(data <10000);
[h x]=hist(data(ind),100);
%p1=plot(x,h);%
p1=plot(x,log(h));
set(p1,'Color','blue','LineWidth',1.0);
xlabel('ADU''s');
ylabel('frequency');
title('data histogram');


end

function data = init_crop(data,nnc) 

%% initial cropping

if numel(nnc) == 6
    
    
    %if sum(nnc(1:2)) > 0, data=crop_dim(data,nnc(1:2),1);end
    %if sum(nnc(3:4)) > 0, data=crop_dim(data,nnc(3:4),2);end
    %if sum(nnc(5:6)) > 0, data=crop_dim(data,nnc(5:6),3);end
    
    sz=size(data);
    xs=sz(2);
    ys=sz(1);
    zs=sz(3);
    
    if sum(nnc(1:2)) < 0,
        disp('doing intial x cropping....')
        nncc=0;
        nncc=[abs(nnc(1)),xs-abs(nnc(2))];
        disp(['xs [',num2str(xs),'] --> [',num2str(nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,1);end
    
    if sum(nnc(3:4)) < 0, 
        disp('doing intial y cropping....')
        nncc=0;
        nncc=[abs(nnc(3)),ys-abs(nnc(4))];
        disp(['ys [',num2str(ys),'] --> [',num2str(nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,2);end
   
    
    if sum(nnc(5:6)) < 0, 
        disp('doing intial z cropping....')
        nncc=0;
        nncc=[abs(nnc(5)),zs-abs(nnc(6))];
        disp(['zs [',num2str(zs),'] --> [',num2str(nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,3);end
else
    disp('if initial cropping is required')
    disp('set nnc=[-x0,-x1,-y0,-y1,-z0,-z1], where nnc gives pixels to crop')
    disp('off each dimension otherwise set nnc=[0]')
end

end

function data = remove_aliens(data,aliens)
%% determine aliens
a_c=0;      %alien counter for the loop
disp(' ')
disp('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>')
disp('determining alien (<>) removal....')

nx=size(data);

if sum(sum(aliens)) == 0,
    disp('no alien removal....')
else
    disp('aliens detected, preparing to remove....')
    tri=numel(aliens)/3.0;
    nal=numel(aliens)/6.0;          %6 comes from needing 2 pairs of 3 coords
    
    if mod(nal,1) == 0,             %do removal if there are the correct number of paramters
        xyz=reshape(aliens,3,nal*2);        %puts it into the format  [x0,x1,x2...]
                                            %                         [y0,y1,y2...]
                                            %                         [z0,z1,z2...]    
        for hh = 1:2:tri,               %do every second to get the pairs of points
            a_c=a_c+1;
            disp(['removing alien - [',num2str(a_c),'/',num2str(nal),']'])
            
            xx=sort(xyz(1,hh:hh+1));  %x points, automatcially gets the corect order
            yy=sort(xyz(2,hh:hh+1));  %y points  low to high
            zz=sort(xyz(3,hh:hh+1));  %z points
            
            if min(xx) == -1,xx=[max(xx),nx(2)];end
            if min(yy) == -1,yy=[max(yy),nx(1)];end
            if min(zz) == -1,zz=[max(zz),nx(3)];end
            
            
            disp(['[xstart xfinish] = [',num2str(xx),']'])
            disp(['[ystart yfinish] = [',num2str(yy),']'])
            disp(['[zstart zfinish] = [',num2str(zz),']'])
            
            data(yy(1):yy(2),xx(1):xx(2),zz(1):zz(2))=0;
            
        end                                    
    else
        disp(' ')
        disp('#<># ERROR #<># ERROR #<># ERROR #<># ERROR #<># ERROR #<>#')
        disp(' ')
        disp('incorrect number of [x,y,z] pairs.  require 6 points per alien.')
        disp('no alien removal today....')
        disp(' ')
        disp('#<># ERROR #<># ERROR #<># ERROR #<># ERROR #<># ERROR #<>#')
    end
end
disp(' ')
disp('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>') 

end

function data = init_pad(data,nnc)

if numel(nnc) == 6
    
    if sum(nnc(1:2)) > 0, 
        disp('doing intial x padding....')
        data=padarray(data,[0 abs(nnc(1)) 0],0,'pre');
        data=padarray(data,[0 abs(nnc(2)) 0],0,'post');
    end
    
    if sum(nnc(3:4)) > 0, 
        disp('doing intial y padding....')
        data=padarray(data,[abs(nnc(3)) 0 0],0,'pre');
        data=padarray(data,[abs(nnc(4)) 0 0],0,'post');
    end
    
    if sum(nnc(5:6)) > 0, 
        disp('doing intial z padding....')
        data=padarray(data,[0 0 abs(nnc(5))],0,'pre');
        data=padarray(data,[0 0 abs(nnc(6))],0,'post');
    end
    
else
    disp('if initial pading is required')
    disp('set nnc=[+x0,+x1,+y0,+y1,+z0,+z1], where nnc gives pixels to pad')
    disp('onto each dimension otherwise set nnc=[0]')
end


end

function data = fft_pad(data,bin)

%% set an array of values that it will automtically pad to to give 
%a good number for the FFT, i.e powers of 2 or PX2^N, where P is a low
%numbered prime i.e <7 or 9
%pref=[32,48,64,80,96,128,144,160,192,256,320,512];
pow_2=2.^(1:10);
pref=sort([2.^(6:10),2.^(5:10)*3,2.^(4:10)*5,2.^(3:10)*9]);%sort([2.^(6:10),pow_2*3,pow_2*5,pow_2*9]);

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

    %disp(['Array size after binning [x,y,z] - [',num2str([sx,sy,sz]),']'] )

    while sum(sx == bin(1)*pref) == 0,sx=sx+1;end
    while sum(sy == bin(2)*pref) == 0,sy=sy+1;end
    while sum(sz == bin(3)*pref) == 0,sz=sz+1;end

    nn=[sx,sy,sz];
    disp(['Padding for FFT [x,y,z] - [',num2str([sx,sy,sz]),']'])
    disp(['After binning [x,y,z] - [',num2str([sx/bin(1),sy/bin(2),sz/bin(3)]),']'])

    data=zero_pad_ver2(data,nn(1),nn(2),nn(3));
end
if numel(dsz) == 2
    sx=dsz(2);
    sy=dsz(1);

    %disp(['Array size after binning [x,y] - [',num2str([sx,sy]),']'] )

    while sum(sx == bin(1)*pref) == 0,sx=sx+1;end
    while sum(sy == bin(2)*pref) == 0,sy=sy+1;end

    nn=[sx,sy];
    disp(['Padding for FFT [x,y] - [',num2str([sx,sy]),']'])
    disp(['After binning [x,y,z] - [',num2str([sx/bin(1),sy/bin(2)]),']'])


    data=zero_pad_ver2(data,nn(1),nn(2));
end

disp('Complete.....')
disp(' ')

%disp('----------------------------------------------------------------')


end

function saturation_warning(file)

disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('!! WARNING !! WARNING !! WARNING !!')
disp('!! Multiple maximum values found !!')
disp('!!   possibly saturated data     !!')
disp('***********************************')  
disp(['** check ',file,' **'])
disp('***********************************')  



end

