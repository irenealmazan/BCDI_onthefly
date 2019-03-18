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
% the last file of your .spe data set can also be used as a bg (for example
% if no bg's were taken and the last file is essentially a bg), just set
% back={'last'}.  you can provide some bg's and use the last eg.
% files={'file1.spe','file2.spe'}
% back={'bg.spe','last'}
% will use the bg for the first file and the last frame of the second file
% data set as the bg.
%
%
%for a series of .tif files, it is assumed these have come from a pilatus and the 3d
%data set is made up of n 2d images.  no bg 
%subtraction is done.  just enter the first and last file names and the
%number of accuulations per theta position (file_params.nfiles=#) and it
%will make the 3d data set.
%
%if only 1 .tif file is provided it assumes that it is 3d.
%
%alternatively provide a .mat file and it will directly read in a matlab
%file, assumed to be already in 3d.
%
%providing .spe,.tif or.mat will still use mindata do remove values below
%this and will center the data onto the maximum and remove aliens.

%file_params contains other stuff

%file_params.subtract_dc == 1 will subtract the mindata value off the data
%                           rather than just setting everything below this to 0.

%file_params.schot_th = 200 will do a secondary threshold after binning
%                   useful since some schot noise won't get eliminated by the first threshold
%                   but a second higer one will

%file_params.no_center = 1, won't center the data, useful for ptycho.

%file_params.no_fft_pad=1, won't do any fft padding
%file_params.pad_ptych=1; will do fft padding inx and y (i.e for ptycho)
%file_params.no_hist=1; won't plot the histograms
%file_params.bg_mult=1; multiplication to apply to the bg file (if exp time is different)
%file_params.do_2D=1; will return a 2d slice (central)
%file_params.nnc=[-10,-10,0,0,2,4] does initial cropping and padding.  -ve
%                 numbers are crop, +ve are pad.  it goes
%                 nnc=[left,right,top,bottom,start,finish]
%file_params.return_orig_size=1 ; returns the data array with the original
%               size.  default is 0 (don't enforce it to return original size)
%               used for .spe files only at this stage (check though)
%
%file_params.data_shift=[x,y,z] ; will shift the data  default is 0 (not
%                                   supported yet)
%
% check the function bin_crop_center_defaults() which is at the end of this
% file for every option.  some might not have made it to these notes.
%% set defaults if not already
disp('----------------------------------------------------------------')

try
    file_params;
catch
    file_params=[];
end

file_params=bin_crop_center_defaults(file_params);

%%

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
        tmp=char(bgs(1));
        if strcmp(tmp(end-3:end),'last')
            disp('using last frame as background....')
            back0=data0(:,:,end);
            
        else
            back0=speread(char(bgs(1)));          %get background image
            disp('Loading background image - [1] ')
            back0=rot3d(back0.data);

            if ndims(back0) == 3,back0=back0(:,:,2);end
            
    
        end
    end   
        
    if numel(nx) == 2,nx=[nx,1];end

    data=data0-data0;
    for qq = 1:nx(3), data(:,:,qq)=flipud(data0(:,:,qq)-file_params.bg_mult*back0);end         %remove bg
    data0=0;
    back0=0;

    ind=(data < mindata);             %keep everything above mindata
    data(ind)=0;

    
    if file_params.subtract_dc == 1
        disp('subtracting mindata from data....')
        data=data-mindata;
        ind=(data < 0);
        data(ind)=0;
        ind=-1;
    end
    
    %pad it first jnc nov 2011
    %data = init_pad(data,nnc);
    
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
                tmp=char(bgs(ww));
                if strcmp(tmp(end-3:end),'last')
                    
                    disp('using last frame as background....')
                    back1=data1(:,:,end);
            
                else

                    disp(['Loading background image - [',num2str(ww),']'])
                    back1=speread(char(bgs(ww)));          %get background image
                    back1=rot3d(back1.data);

                    %if exist('back_no') == 0,back_no=ndims(back1)-1;end
                    %if back_no == 2,back1=back1(:,:,back_no);end

                    if ndims(back1) == 3,back1=back1(:,:,2);end
                    %back1=back1(:,:,2);
                end
            end
            
            for qq = 1:nx(3), data1(:,:,qq)=flipud(data1(:,:,qq)-file_params.bg_mult*back1);end  
           
            ind=( data1 < mindata );
            data1(ind)=0;
            ind=-1;
            %data1=center_array(data1);
            
            if file_params.subtract_dc == 1
               disp('subtracting mindata from data....')
               data1=data1-mindata;
               ind=(data1 < 0);
               data1(ind)=0;
               ind=-1;
            end
            
            %pad it first jnc nov 2011
            %data1 = init_pad(data1,nnc);
            
            [v1 v2 v3]=ind2sub(size(data1),find(data1 == max(max(max(data1)))));  %get center of next data set
            
            if numel(v1) ~= 1,
                saturation_warning(char(files(ww)));
                v1=v1(1);
                v2=v2(1);
                v3=v3(1);
            end
            
            if u3-v3 ~= 0,
                disp('Data sets don''t'' share common Z center....')
                disp(['Centers off by ',num2str(u3-v3),' pixels....'])
                disp('Zeroing start/end frame to prevent trouble....')
                
                if u3-v3 > 0,data1(:,:,end-abs(u3-v3)+1:end)=0;end
                if u3-v3 < 0,data1(:,:,1:abs(u3-v3))=0;end
                
            end
                
            data1=circshift(data1,[u1-v1,u2-v2,u3-v3]);
            
            data=data+data1;
            data1=-1;
        end

    end
    
    data=remove_aliens(data,aliens);
    s0=size(data);
    if file_params.pad_ptych ~= 1,disp(['Data raw size [x,y,z] = [',num2str(s0(2)),',',num2str(s0(1)),',',num2str(s0(3)),']']),end
    data = init_pad(data,nnc);
    data = init_crop(data,nnc);
    s01=size(data);
    
    if sum(abs(s0-s01)) ~=0
        disp(['Data initial pad/crop size [x,y,z] = [',num2str(s01(2)),',',num2str(s01(1)),',',num2str(s01(3)),']'])
    end
    
    if file_params.no_fft_pad == 0
        if file_params.pad_ptych == 1
            
            data=fft_pad_ptych(data,[bin,1]);
        else
            data=fft_pad(data,[bin,1]);
        end
        
    else
        
        disp('No FFT padding....')
    
    end
    
    if file_params.no_center == 0,
        data=center_array(data);
        disp('Centering data....')
    else
        disp('Not centering data....')
    end

    if file_params.photon_heal == 1,
        disp('Photon healing....')
        data=photon_healing(data,ph_th,ph_sig);
    end
    
    if file_params.return_orig_size == 1
        if sum(abs(s0-size(data))) ~= 0
           disp('Returning data with original size....')
           
           s02=size(data);
           ncp=s0-s02;
           nnc1=[floor(ncp(2)/2),ceil(ncp(2)/2),floor(ncp(1)/2),ceil(ncp(1)/2),floor(ncp(3)/2),ceil(ncp(3)/2)];
           data = init_pad(data,nnc1);
           data = init_crop(data,nnc1);
        end
    end
    
    
end

if isempty(strfind(lower(char(files{1})),'.tif')) == 0
    
    temp = imfinfo(char(files{1}), 'tif');
    
    if numel(temp) == 1
        
        disp('Loading multiple .tif data files to stack....')
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
       
        if file_params.no_fft_pad == 0
            if file_params.pad_ptych == 1

                data=fft_pad_ptych(data,[bin,1]);
            else
                data=fft_pad(data,[bin,1]);
            end
        
        else
            disp('No FFT padding....')
        end
        
        if file_params.no_center ==0,
            data=center_array(data);
            disp('Centering data....')
        else
            disp('Not centering data....')
        end
        
        ind=( data < mindata );
        data(ind)=0;
    
    else
        
        nfiles=size(files);             %get the number of files
        nfiles=max(nfiles);
        
        data_s=0;                       %will have the summed data (mult data sets)
        
        disp(['Loading [',num2str(nfiles),'] .tif stacks....'])
        
        for ww = 1:nfiles
            
            fnames=char(files{ww});
            disp('')
            disp('          ############################')
            disp(['          # Loading .tif stack [',num2str(ww),'/',num2str(nfiles),'] #'])
            disp('          ############################')
            
            temp = imfinfo(char(files{ww}), 'tif');
            
            num_images=numel(temp);
            for qq = 1:num_images

                ff0 = imread(fnames,qq);

                if qq == 1, data=zeros([size(ff0),num_images]);end

                data(:,:,qq) = ff0;

            end

            data=remove_aliens(data,aliens);
            data = init_pad(data,nnc);
            data = init_crop(data,nnc);
            
            if file_params.no_fft_pad == 0
                if file_params.pad_ptych == 1

                    data=fft_pad_ptych(data,[bin,1]);
                else
                    data=fft_pad(data,[bin,1]);
                end

            else
                disp('No FFT padding....')
            end
            if file_params.no_center ==0,
                data=center_array(data);
                disp('Centering data....')
            else
                disp('Not centering data....')
            end
            ind=( data < mindata );
            data(ind)=0;
            data_s=data_s+data;                 %summ data sets
        end
        disp('Done....')
        data=data_s;
        data_s=0;
    end

end

if isempty(strfind(lower(char(files{1})),'.mat')) == 0
    
    disp('Loading MATLAB data file....')
    
    load(char(files{1}));
    data_a=data;                %assumes data is called 'data'
    
    if numel(bgs) == 0,back0=0;else   %check if a bg files is supplied
        tmp=char(bgs(1));
        disp('Loading MATLAB background file....')
        load(tmp);
        back0=data;
        data=0;
        if ndims(back0) == 3, back0=back0(:,:,2);end
        
        if ndims(data_a) == 3,nz=size(data_a,3);else nz=1;end
        
        if nz > 1
           for qq = 1:nz,
               data_a(:,:,qq)=data_a(:,:,qq)-file_params.bg_mult*back0;
           end
           data=data_a;
           data_a=0;
        end
        
        if nz == 1
           data=data_a-file_params.bg_mult*back0; 
        end
        
    end
    
    nx=size(data);
    disp(['Values below -[',num2str(mindata),']- will be set to 0...'])
    disp(' ')
    ind=( data < mindata );
    data(ind)=0;
    
    if file_params.subtract_dc == 1
       disp('subtracting mindata from data....')
       data=data-mindata;
       ind=(data < 0);
       data(ind)=0;
    end
   
    data=remove_aliens(data,aliens);
    data = init_pad(data,nnc);
    data = init_crop(data,nnc);
    
    if file_params.no_fft_pad == 0
        if file_params.pad_ptych == 1

            data=fft_pad_ptych(data,[bin,1]);
        else
            data=fft_pad(data,[bin,1]);
        end

    else
        disp('No FFT padding....')
    end
    if file_params.no_center==0,
        data=center_array(data);
        disp('Centering data....')
    else
        disp('Not centering data....')
    end
    
    
end

if isempty(strfind(lower(char(files{1})),'.h5')) == 0

    
    nfiles=size(files);             %get the number of files
    nfiles=max(nfiles);

    data_s=0;                       %will have the summed data (mult data sets)

    disp(['Loading [',num2str(nfiles),'] .h5 stacks....'])

    for ww = 1:nfiles

        fnames=char(files{ww});
        disp('')
        disp('          ############################')
        disp(['          # Loading .h5 stack [',num2str(ww),'/',num2str(nfiles),'] #'])
        disp('          ############################')

        %read the data
        data=hdf5read(fnames,'data');

        %read the background
        if numel(bgs) == 0,back0=0;else   %check if a bg files is supplied
            tmp=char(bgs(1));
            disp('Loading .h5 background file....')
            back0=hdf5read(tmp,'data');             %read in bg
         
            if ndims(back0) == 3, back0=back0(:,:,2);end  %check if mult bg's
            if ndims(data) == 3,nz=size(data,3);else nz=1;end   %check if 2d

            if nz > 1                   %if mult bg's tkae 2nd one
               for qq = 1:nz,
                   data(:,:,qq)=data(:,:,qq)-file_params.bg_mult*back0;
               end
            end

            if nz == 1
               data=data-file_params.bg_mult*back0; 
            end
        end
        
        data=remove_aliens(data,aliens);
        data = init_pad(data,nnc);
        data = init_crop(data,nnc);

        if file_params.no_fft_pad == 0
            if file_params.pad_ptych == 1

                data=fft_pad_ptych(data,[bin,1]);
            else
                data=fft_pad(data,[bin,1]);
            end

        else
            disp('No FFT padding....')
        end
        if file_params.no_center ==0,
            data=center_array(data);
            disp('Centering data....')
        else
            disp('Not centering data....')
        end
        ind=( data < mindata );
        data(ind)=0;
        data_s=data_s+data;                 %summ data sets
    end
    disp('Done....')
    data=data_s;
    data_s=0;


end

%%
nx=size(data);


if file_params.no_hist == 0
    plot_hist_data(data,0)
end

%pad the array so that it bins exactly
x0=nx(2);
y0=nx(1);

if max(size(nx)) == 2,nx=[nx,1];end

disp(' ')
disp('Resizing data....')
disp(['Current data size [x,y,z] - [',num2str([x0,y0,nx(3)]),']'])

while mod(x0,bin(1)) ~=0,x0=x0+1;end
while mod(y0,bin(2)) ~=0,y0=y0+1;end


data=padarray(data,[(y0-nx(1)),(x0-nx(2)),0],0,'pre');

%data=imresize(center_array(data),[y0/bin(2),x0/bin(1)],'nearest');
data_new=zeros([y0/bin(2),x0/bin(1),nx(3)]);

for qq = 1:nx(3),     
    if bin(1) == 1 && bin(2) == 1
        data_new=data;
        if qq == 1,disp('No binning since bx = 1 and by = 1'),end
    else
        data_new(:,:,qq)=box_interp(data(:,:,qq),bin(1),bin(2));
    end
end

data=data_new;
data_new=0;

disp(['Array size after binning [x,y,z] - [',num2str([x0/bin(1),y0/bin(2),nx(3)]),']'])

try
    disp(['doing secondary thresholding - [',num2str(file_params.schot_th),']'])
    ind=(data < file_params.schot_th);
    data(ind)=0;
end

%adjust orientation of necessary
data=adjust_orientation(data,file_params.det_orient);


if file_params.do_2D == 1
    
    if ndims(data) > 2
       data=data(:,:,round(size(data,3)/2));
    end
    disp('Returning 2D (central) slice....') 
end

if file_params.no_hist == 0;
    plot_hist_data(data,1)
end


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
    if ndims(data) ~= 2,zs=sz(3);end
    
    if sum(nnc(1:2)) < 0,
        disp('doing intial x cropping....')
        nncc=0;
        nncc=[abs(nnc(1))+1,xs-abs(nnc(2))];
        disp(['xs [',num2str(xs),'] --> [',num2str(1+nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,1);end
    
    if sum(nnc(3:4)) < 0, 
        disp('doing intial y cropping....')
        nncc=0;
        nncc=[abs(nnc(3))+1,ys-abs(nnc(4))];
        disp(['ys [',num2str(ys),'] --> [',num2str(1+nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,2);end
   
    
    if sum(nnc(5:6)) < 0, 
        disp('doing intial z cropping....')
        nncc=0;
        nncc=[abs(nnc(5))+1,zs-abs(nnc(6))];
        disp(['zs [',num2str(zs),'] --> [',num2str(1+nncc(2)-nncc(1)),']'] )
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

sz=size(data);
xs=sz(2);
ys=sz(1);

try
 zs=sz(3);
end
 
if numel(nnc) == 6
    
    if sum(nnc(1:2)) > 0, 
        disp('doing intial x padding....')
        data=padarray(data,[0 abs(nnc(1)) 0],0,'pre');
        data=padarray(data,[0 abs(nnc(2)) 0],0,'post');
        disp(['xs [',num2str(xs),'] --> [',num2str(xs+nnc(2)+nnc(1)),']'] )
    end
    
    if sum(nnc(3:4)) > 0, 
        disp('doing intial y padding....')
        data=padarray(data,[abs(nnc(3)) 0 0],0,'pre');
        data=padarray(data,[abs(nnc(4)) 0 0],0,'post');
        disp(['ys [',num2str(ys),'] --> [',num2str(xs+nnc(4)+nnc(3)),']'] )
    end
    
    if sum(nnc(5:6)) > 0, 
        disp('doing intial z padding....')
        data=padarray(data,[0 0 abs(nnc(5))],0,'pre');
        data=padarray(data,[0 0 abs(nnc(6))],0,'post');
        disp(['zs [',num2str(zs),'] --> [',num2str(zs+nnc(5)+nnc(6)),']'] )
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
pref=sort([2.^(6:16),2.^(5:10)*3,2.^(4:10)*5,2.^(3:10)*9]);%sort([2.^(6:10),pow_2*3,pow_2*5,pow_2*9]);

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

    data=zero_pad_ver2(data,nn(1),nn(2),nn(3));
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

function data = fft_pad_ptych(data,bin)

%% set an array of values that it will automtically pad to to give 
%a good number for the FFT, i.e powers of 2 or PX2^N, where P is a low
%numbered prime i.e <7 or 9
%pref=[32,48,64,80,96,128,144,160,192,256,320,512];
pow_2=2.^(1:10);
pref=sort([2.^(6:16),2.^(5:10)*3,2.^(4:10)*5,2.^(3:10)*9]);%sort([2.^(6:10),pow_2*3,pow_2*5,pow_2*9]);

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

    while sum(sx == bin(1)*pref) == 0,sx=sx+1;end
    while sum(sy == bin(2)*pref) == 0,sy=sy+1;end
    %while sum(sz == bin(3)*pref) == 0,sz=sz+1;end

    nn=[sx,sy,sz];
    disp(['Padding for FFT [x,y,z] - [',num2str([sx,sy,sz]),']'])
    disp(['After binning it will be, [x,y,z] - [',num2str([sx/bin(1),sy/bin(2),sz/bin(3)]),']'])

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

function healed=photon_healing(data,th,sig)


auto=ifftshift(ifftn(fftshift(data)));

%a_sup=shrink_wrap(abs(auto).^.5,th,sig);
sx=100;
nn=size(auto);



healed=real(ifftshift(fftn(fftshift(auto))));

end

function file_params=bin_crop_center_defaults(file_params)

try
    file_params.no_center;
catch
    file_params.no_center=0; %default to centering (0=yes, 1=no)
end

try
    file_params.no_fft_pad;
catch
    file_params.no_fft_pad=0; %default to fft padding (0=yes, 1=no)
end

try
    file_params.pad_ptych;
catch
    file_params.pad_ptych=0; %default is to do 3d padding rather than 2d for ptych
end

try
    file_params.no_hist;
catch
    file_params.no_hist=0;
end

try
    file_params.bg_mult;
catch
    file_params.bg_mult=1;
end
try
    file_params.do_2D;
catch
    file_params.do_2D=0;
end
try
    file_params.subtract_dc;
catch
    file_params.subtract_dc=0;
end

try
    file_params.photon_heal;
    ph_th=.05;
    ph_sig=5;
catch
    file_params.photon_heal=0;
end

try
    file_params.det_orient;
catch
    file_params.det_orient=[];
end
try
    file_params.return_orig_size;
catch
    file_params.return_orig_size=0;
end
try
    file_params.data_shift
catch
    file_params.data_shift=0;
end

end

function data=adjust_orientation(data,orient);


if isempty(orient) == 0
   disp(' ') 
   disp('/|\_/|\_/|\_/|\_/|\_/|\_/|\_/|\_/|\_/|') 
   switch orient
       
       case {'LCLS','lcls'}
           disp('Orientating data for LCLS setup....')
           data=flipdim(flipdim(data,1),2);
           
       case {[],'none','NONE'}
           disp('Leaving unchanged....')
           
       otherwise
           disp('Unknown orientation, leaving unchanged....')
   end
   
   disp('/|\_/|\_/|\_/|\_/|\_/|\_/|\_/|\_/|\_/|') 
   disp(' ') 
    
else
   
    disp('No change to data orientation....')
    
end

end

function plot_hist_data(data,bef_aft)


if bef_aft == 1,
    disp(' ')
    figure(90)
    disp('Plotting histogram after binning....'),
    ind=(data <10000);
    [h x]=hist(data(ind),100);
    %p1=plot(x,h);%
    p1=plot(x,log(h));
    set(p1,'Color','blue','LineWidth',1.0);
    xlabel('ADU''s');
    ylabel('frequency');
    title('data histogram');
else
    disp(' ')
    figure(89)
    disp('Plotting histogram before binning....')
    ind=(data <2000);
    [h x]=hist(data(ind),100);
    %p1=plot(x,h);%
    p1=plot(x,log(h));
    ymin=log(h(2:end));
    ymin=min(ymin(ymin>0));
    %axis([0 2000 ymin max(log(h(2:end)))])
    set(p1,'Color','blue','LineWidth',1.0);
    xlabel('ADU''s');
    ylabel('frequency');
    title('data histogram');
end

end