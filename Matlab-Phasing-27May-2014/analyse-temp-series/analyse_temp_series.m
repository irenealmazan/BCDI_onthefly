function  analyse_temp_series(top_dir,prefix,maxy)
%jclark
% to analyse a temperature series of CDI
% it is assumed the reconstruction to be analysed will be an average
% from several random starts.  with the summation performed using
% sum_align_reconstructions
% top_dir should be the dir that has the temperatures below it
% e.g. top_dir/T030/rand-starts/CVL/

no_data=0;    %don't do the data, its takes time, use data = 0
avgr=3;       %region to average over, nxn
offset=1;     %plot the lineouts with an offset?, yes=1

shell_hist=0;   %do the shell histograms, yes=1
amp_norm=0;     %normalise amplitudes of reconstruction to that of the first
                %such that max(first)=1, and then
                %sum(second)=sum(first).....
plot_cold=1;    %set =1 to plot the cooled down temperature 30 if present
     
cold_string='a/r';  %s string that is in the file name of the cold temps that distniguishes them from the rest

lab=1;          %load the transformed reconstruction

vtk_out=1;      %output vtk file (will be correctly orientated)

auto_size_plot=0;   %auto size (=0), use all =1

conj_first=0;     %conjugate the reference? (=1 yes, 0=no)

rem_ramp=1;     %remove ramp in sample plane, =1 yes
upsize=3;       %upsample factor

z_phase=1;      %set mean to zero (1=yes)

apply_sup=1;    %multiply by support (vtk output)

hist_range=.5;     %fraction of 2 pi
ph_bins=100;

ph_range=[-pi,pi];  %range for output image

offset_amp=.65;          %offset to plot amp histogam as fraction of max
offset_ph=.6;          %offset to plot ph histogam as fraction of max

vsig=0.5;               %calc values for the volume
vth=0.1;

stuff.bold_plot=0;
stuff.dots_only=0;

sort_numerically=1;%1;     %will sort according to number.  do not use if heating series has coold down temps, will lose the order
nnumbs_out=[2];%[2];       %number of digits in the numbers for save name output, use [] to turn off

string_prefix='EV';%'T';

custom_xaxis=[];%[10,20,50,80,110,140,160,200,230,290,-40];        %use a custom axis, leave as [] for nothing

cxlabel='Time (x20 minutes)';%'Delay time (ps)';%;'Temperature (C)';%'

points=[30,47,36];  %x,y,z, J9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% @
%output the above params
file=[return_analyse_temp_series_script_path(),'analyse_temp_series.m'];
[ str_n ] = load_script_as_text(file,'@');
%%
try 
    prefix;
catch
    prefix='Rec-5x7-*-GAHIO-AVGh-mC-10-5-200-CVl-SW';
end

try
    maxy;
    params.maxy=maxy;
catch
    params.maxy=1000;
end
maxy=params.maxy;

try
    top_dir;
catch
    name_of_this_file='analyse_temp_series';
    dir_file=which(name_of_this_file);    %finds the current folder for the phasing
    top_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
end
try
    nnumbs_out;
catch
    nnumbs_out=[];
end

try
    cold_string;
catch
    cold_string='a/R';
end
%%
% get names of the summed and averaged reconstructions and a paramter file
% for each reconstruction

%string='rand-starts/CVL/CVL-AMP.rec';

if lab ~= 1,string_amp=[prefix,'**-AMP.rec'];else string_amp=[prefix,'**LAB.rec'];end
amp_fs=rdir([top_dir,'**/*/',string_amp]);
n_temps=size(amp_fs,1);                         %get the number of temperatures

string_ph=[prefix,'**-PH.rec'];
ph_fs=rdir([top_dir,'**/*/',string_ph]);

%amp_fs_c=rdir([top_dir,string_prefix,'*a*/*/',string_amp]); %get the names of the cold ones, if any
amp_fs_c=rdir([top_dir,string_prefix,'*',cold_string,'*/*/',string_amp]); %get the names of the cold ones, if any
n_temps_cold=size(amp_fs_c,1);                 %get the number of cold temps, an use this as a sanity check later on
check=zeros([1,n_temps]);

string_params=[prefix,'**PARAMS.mat'];
pms_fs=rdir([top_dir,'**/*/',string_params]);

for qq=1:n_temps
    if numel(regexp(char(amp_fs(qq).name),'a/r','match')) > 0,check(qq)=1;end
    %if numel(regexp(char(amp_fs(qq).name),'a/r','match')) >
    %0,check(qq)=1;end  careful of the lower or upper case a/r or a/R
end

if sum(check) ~= n_temps_cold,disp('ERROR determining cold temperatures....'), else disp('Cold temperatures determined....'),end

%%

warm_ind=find(check == 0);
for qq=1:numel(warm_ind)

    temp_a_fs(qq).name=char(amp_fs(warm_ind(qq)).name);
    temp_ph_fs(qq).name=char(ph_fs(warm_ind(qq)).name);
    temp_pm_fs(qq).name=char(pms_fs(warm_ind(qq)).name);
end

if plot_cold == 1, %=0 -> don't plot cold
    
    cold_ind=reverse(find(check == 1));
    for qq=1:numel(cold_ind)

        temp_a_fs(qq+numel(warm_ind)).name=char(amp_fs(cold_ind(qq)).name);
        temp_ph_fs(qq+numel(warm_ind)).name=char(ph_fs(cold_ind(qq)).name);
        temp_pm_fs(qq+numel(warm_ind)).name=char(pms_fs(cold_ind(qq)).name);
        
    end
end

amp_fs=temp_a_fs;
ph_fs=temp_ph_fs;
pms_fs=temp_pm_fs;
n_temps=max(size(amp_fs));                         %get the number of temperatures (without cold)

%sort the numbers, be careful if you want cooling doen temps in order, =0
%for this
if sort_numerically == 1
   disp(' ')
   disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
   disp('SORTING SERIES....BE CAREFUL IF PLOTTING COOLING DOWN AS WELL....')
   disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
   disp(' ')
   if isempty(custom_xaxis) == 1
       for rr=1:n_temps,numbs_init(rr)=extract_number_from_string(char(amp_fs(rr).name),[],string_prefix);end
       [BB IIX]=sort(numbs_init);
   else
       [BB IIX]=sort(custom_xaxis);
       custom_xaxis=custom_xaxis(IIX);
   end
   amp_fs=amp_fs(IIX);
   ph_fs=ph_fs(IIX);
   pms_fs=pms_fs(IIX);

end

%%
% load the first one to get the dimensions of the reconstriuction
if lab ~=1,first = load_rec(strtrim(char(amp_fs(1).name)),strtrim(char(ph_fs(1).name)));
else
    load(strtrim(char(amp_fs(1).name)),'-mat')
    first=array;
    array=[];
end

if conj_first == 1,
    disp('Conjugating and reflecting first....')
    first=conj_reflect(first);end


if rem_ramp == 1,first = remove_ramp_pn_ups(first,upsize);end
first=first/max(abs(first(:)));

if z_phase == 1,first=zero_phase(first,0);end

tot_first=sum(abs(first(:)));
load(pms_fs(1).name);

if no_data == 0, data=load_rec_data(params,[],top_dir);end

% make an array to store the all
temp_series=zeros([size(first),n_temps]);
vol_series=zeros([size(first),n_temps]);

temp_series(:,:,:,1)=first;
ref=abs(first);

%some prelims
tot_n=numel(abs(first));
th=0.2;
sig=0.5;
volume=0;
sphericity=0;
nxcc=0;
strain=0;
strain_sig=0;
std_dev=0;
numbers=0;
skew_ph=0;
numbers_string={};
params.colors=['r','b'];        %used for the histogram colors
nshells=2;                      %number of shells for the histogram

save_dir=[top_dir,'T-analysis'];
if isdir(save_dir) == 0,mkdir(save_dir);end
vtk_dir=[save_dir,'/vtks'];
if isdir(vtk_dir) == 0,mkdir(vtk_dir);end

xy_dir=[save_dir,'/xy'];
if isdir(xy_dir) == 0,mkdir(xy_dir);end
xz_dir=[save_dir,'/xz'];
if isdir(xz_dir) == 0,mkdir(xz_dir);end
zy_dir=[save_dir,'/zy'];
if isdir(zy_dir) == 0,mkdir(zy_dir);end


%get temperature from directory name or use custom xaxis
if isempty(custom_xaxis) ==  1
    nums=char(regexp(char(amp_fs(1).name),[string_prefix,'\d+'],'match'));  %nums=char(regexp(char(amp_fs(1).name),'T\d+','match'));
    nums=nums(numel(string_prefix)+1:end);                                  %nums=nums(2:end);
else
    nums=num2str(custom_xaxis(1)); 
    disp('Using custom xaxis....')
end
    
if isempty(nnumbs_out) ~= 1,nums = check_numlength(nums,nnumbs_out);end

numbers_string(1)={nums};
numbers(1)=str2num(nums);

%save the names of the files, cold and hot
output_file_names([save_dir,'/File_names.txt'],amp_fs)
%save the paramters
output_string([save_dir,'/temp-series-params.txt'],str_n)

%calc the first volume and cor-coeff and amp values spread (std dev)
vol_sup=shrink_wrap(abs(first),vth,vsig);       %volume over which to calc stuff
vol_series(:,:,:,1)=vol_sup;
vol_ref=numel(find(vol_sup > 0));
volume(1)= numel(find(vol_sup > 0))/vol_ref;   %as % of total
nxcc(1)=1.0;%max(max(max(normxcorr3(abs(first), abs(first), 'same')))); %normalised cross-corr
stdd=abs(first(find(vol_sup > 0) ));
std_dev(1)=std(stdd(:));
amp_mean(1)=mean(stdd(:));
sphericity(1)=calc_sphericity(vol_sup);

%output a vtk into the dir if specified
if vtk_out == 1
    vtkname=[vtk_dir,'/',string_prefix,nums,'-Amp-Phase.vtk'];
    
    if apply_sup == 1,vtkarray=flip_all_dim(first.*vol_sup);
    else vtkarray=flip_all_dim(first);end
    
    savevtk2scalar(abs(vtkarray),vtkname,angle(vtkarray),1)
    vtkarray=[];
end

%calc corr coeff for data
if no_data == 0, dxcc(1)=max(max(max(normxcorr3(data, data, 'same'))));end %normalised cross-corr data

% calc the 'strain'.  it is the integrated phase, would need to be
% unwrapped for large phase
ph=atan2(imag(first),real(first));
ph=ph-mean( ph( find(ph.*vol_sup ~=0)) );
strain(1)=sum(sum(sum(abs(vol_sup.*ph))))/numel(find(vol_sup > 0));
strain_sig(1)=std( ph(find(ph.*vol_sup ~=0)));
skew_ph(1)=skewness( ph(find(ph.*vol_sup ~=0)));


%output phase and amplitude and volume support
save_dirs={[xy_dir,'/',nums,'I-'],[xz_dir,'/',nums,'I-'],[zy_dir,'/',nums,'I-']};
save_dirsS={[xy_dir,'/',nums,'S-'],[xz_dir,'/',nums,'S-'],[zy_dir,'/',nums,'S-']};
asp = plot_amp_phase(first,[save_dirs],ph_range,auto_size_plot);%.*vol_sup
plot_amp_phase(vol_sup,[save_dirsS],ph_range,asp);
%asp = plot_amp_phase(first,[save_dirs,'/',nums],ph_range,auto_size_plot);%.*vol_sup
%plot_amp_phase(vol_sup,[save_dirs,'/',nums,'S'],ph_range,asp);

%calc histogram for amp + phase
params.prefix='-A-';
params.colors=['r','b'];
params.xaxis='Amplitude Value';
if shell_hist==1;
density_shells_histogram(abs(first),vol_sup,nshells,[save_dir,'/T',char(numbers_string(1)),'/'],params);end

params.prefix='-P-';
params.colors=['c','b'];        %used for the histogram colors
params.xaxis='Phase (2 \pi)';
if shell_hist == 1;
density_shells_histogram((ph+pi)/2/pi,vol_sup,nshells,[save_dir,'/T',char(numbers_string(1)),'/'],params);end

%calc some line-outs
yl=size(first,1)/2-(avgr-1)/2;
yh=size(first,1)/2+(avgr-1)/2;
zl=size(first,3)/2-(avgr-1)/2;
zh=size(first,3)/2+(avgr-1)/2;
xl=size(first,2)/2-(avgr-1)/2;
xh=size(first,2)/2+(avgr-1)/2;
    
%do first temp amps
x_amp_lines(:,1)=mean(mean(abs(first(yl:yh,:,zl:zh)),3),1);
y_amp_lines(:,1)=mean(mean(abs(first(:,xl:xh,zl:zh)),3),2);
z_amp_lines(:,1)=mean(mean(abs(first(yl:yh,xl:xh,:)),2),1);

x_ph_lines(:,1)=mean(mean((ph(yl:yh,:,zl:zh)).*vol_sup(yl:yh,:,zl:zh),3),1);
y_ph_lines(:,1)=mean(mean((ph(:,xl:xh,zl:zh)).*vol_sup(:,xl:xh,zl:zh),3),2);
z_ph_lines(:,1)=mean(mean((ph(yl:yh,xl:xh,:)).*vol_sup(yl:yh,xl:xh,:),2),1);

xls=size(x_amp_lines(:,1),1);
yls=size(y_amp_lines(:,1),1);
zls=size(z_amp_lines(:,1),1);

%track the point
do_point=0;
if exist('points')
   if sum(points(:)) ~= 0
        do_point=1;
        track_point(:,1)=first(points(2),points(1),points(3));
   end
end

for qq = 2:n_temps,
    
    %load the other data and calc cc
    load(pms_fs(qq).name);
    
    if no_data == 0, data_n=load_rec_data(params,[],top_dir);end
    
    if no_data == 0, dxcc(qq)=max(max(max(normxcorr3(data, data_n, 'same'))));end %normalised cross-corr data
    
    
    %load other temp rec and get phase and amp
    if lab ~= 1
        pn=load_rec(strtrim(char(amp_fs(qq).name)),strtrim(char(ph_fs(qq).name)));
    else
        load(strtrim(char(amp_fs(qq).name)),'-mat')
        pn=array;
        array=[];
    end
    %resize in the detector
    pn=ResizeFFt(ref,pn);
    if rem_ramp == 1,pn = remove_ramp_pn_ups(pn,upsize);end
    amp=abs(pn);
    ph=atan2(imag(pn),real(pn));
    
    %align the chirality with the first one
    cnj_rf=is_conj_ref(ref,amp);
            
    if cnj_rf ~=0,crystal=conj_reflect(pn);else
        crystal=(pn);end
    
    crystal=crystal/max(abs(crystal(:)));
    if z_phase == 1,crystal=zero_phase(crystal,0);end
    
    if amp_norm == 1,crystal=crystal*tot_first/sum(abs(crystal(:)));end
    
    amp=abs(crystal);

    disp(max(amp(:)))
    
    %center the reconstructions
    [h k l]=register_3d_reconstruction(ref,amp);

    crystal=(sub_pixel_shift(crystal,h,k,l));
    
    temp_series(:,:,:,qq)=crystal;
    
    %calc volume + nxcc
    vol_sup=shrink_wrap(abs(crystal),vth,vsig);
    vol_series(:,:,:,qq)=vol_sup;
    
    volume(qq)= numel(find(vol_sup > 0))/vol_ref;
    nxcc(qq)=max(max(max(normxcorr3(abs(temp_series(:,:,:,1)), abs(temp_series(:,:,:,qq)), 'same'))));
    sphericity(qq)=calc_sphericity(vol_sup);
    
    %calc some phase stuff, set mean to pi
    ph=atan2(imag(crystal),real(crystal));
    ph=ph-mean( ph(find(ph.*vol_sup ~=0)));
    strain(qq)=sum(sum(sum(abs(vol_sup.*ph))))/numel(find(vol_sup > 0));
    strain_sig(qq)=std( ph(find(ph.*vol_sup ~=0)));
    skew_ph(qq)=skewness( ph(find(ph.*vol_sup ~=0)));
  
    stdd=abs(crystal(find(vol_sup > 0) ));
    std_dev(qq)=std(stdd(:));
    amp_mean(qq)=mean(stdd(:));
    
    if do_point == 1,track_point(:,qq)=crystal(points(2),points(1),points(3));end
    
    %get temp from directory name or use custom xaxis    
    if isempty(custom_xaxis) == 1
        nums=char(regexp(char(amp_fs(qq).name),[string_prefix,'\d+'],'match'));     %nums=char(regexp(char(amp_fs(qq).name),'T\d+','match'));
        nums=nums(numel(string_prefix)+1:end);          %nums=nums(2:end);
    else
        nums=num2str(custom_xaxis(qq));
        disp('Using custom xaxis....')
    end
    
    if isempty(nnumbs_out) ~= 1,nums = check_numlength(nums,nnumbs_out);end
    
    numbers_string(qq)={nums};
    numbers(qq)=str2num(nums);
    
   
    %output amp and phase
    if plot_cold == 1
        if qq > (n_temps-sum(check)),nums_a=[nums,'A'];else nums_a=nums;end
    else nums_a=nums;end
    
    %output a vtk into the dir if specified
    if vtk_out == 1
        vtkname=[vtk_dir,'/',string_prefix,nums_a,'-Amp-Phase.vtk'];
        %vtkarray=flip_all_dim(crystal);
        if apply_sup == 1,vtkarray=flip_all_dim(crystal.*vol_sup);
        else vtkarray=flip_all_dim(crystal);end
        savevtk2scalar(abs(vtkarray),vtkname,angle(vtkarray),1)
        vtkarray=[];
    end

    save_dirs={[xy_dir,'/',nums_a,'I-'],[xz_dir,'/',nums_a,'I-'],[zy_dir,'/',nums_a,'I-']};
    save_dirsS={[xy_dir,'/',nums_a,'S-'],[xz_dir,'/',nums_a,'S-'],[zy_dir,'/',nums_a,'S-']};
    plot_amp_phase(crystal,[save_dirs],ph_range,asp);%.*vol_sup
    plot_amp_phase(vol_sup,[save_dirsS],ph_range,asp);
    %plot_amp_phase(crystal,[save_dirs,'/',nums_a],ph_range,asp);%.*vol_sup
    %plot_amp_phase(vol_sup,[save_dirs,'/',nums_a,'S'],ph_range,asp);
    close all
    
    %histogram of amp values
    params.prefix='-A-';
    params.colors=['r','b'];
    params.xaxis='Amplitude Value';
    if shell_hist == 1
        density_shells_histogram(abs(crystal),vol_sup,nshells,[save_dir,'/T',char(numbers_string(qq)),'/'],params);end

    %hist of phase values
    params.prefix='-P-';
    params.colors=['c','b'];        %used for the histogram colors
    params.xaxis='Phase (2 \pi)';
    if shell_hist == 1
        density_shells_histogram((ph+pi)/2/pi,vol_sup,nshells,[save_dir,'/T',char(numbers_string(qq)),'/'],params);end

    %line-outs
    x_amp_lines(:,qq)=mean(mean(abs(crystal(yl:yh,:,zl:zh)),3),1);
    y_amp_lines(:,qq)=mean(mean(abs(crystal(:,xl:xh,zl:zh)),3),2);
    z_amp_lines(:,qq)=mean(mean(abs(crystal(yl:yh,xl:xh,:)),2),1);

    x_ph_lines(:,qq)=mean(mean((ph(yl:yh,:,zl:zh)).*vol_sup(yl:yh,:,zl:zh),3),1);
    y_ph_lines(:,qq)=mean(mean((ph(:,xl:xh,zl:zh)).*vol_sup(:,xl:xh,zl:zh),3),2);
    z_ph_lines(:,qq)=mean(mean((ph(yl:yh,xl:xh,:)).*vol_sup(yl:yh,xl:xh,:),2),1);

    
    amp=0;
    ph=0;
    pn=0;
    
end

%make 0's in phase lines =NaN so they aren't plotted
x_ph_lines(x_ph_lines == 0)=NaN;
y_ph_lines(y_ph_lines == 0)=NaN;
z_ph_lines(z_ph_lines == 0)=NaN;


% put things we want saved into params structure
params.amp_names=amp_fs;
params.ph_names=ph_fs;
params.volume=volume;
params.sphericity=sphericity;
params.vol_ref=vol_ref;
params.nxcc=nxcc;
params.strain=strain;
params.array_volume=tot_n;
params.array_size=size(crystal);
params.dir=top_dir;


save([save_dir,'/','PARAMS.mat'],'params');

% plot some of the stuff and save images
lw=1.5;
font_size=25;
if plot_cold == 1,stuff.n_temps_cold=n_temps_cold;else stuff.n_temps_cold=0;end

stuff.lw=lw;
stuff.font_size=font_size;
stuff.xlabel=cxlabel;

stuff.ylabel='Volume (\Delta V)';
stuff.save_name=[save_dir,'/Volume-T'];
stuff.color='blue';
plot_generic(numbers,volume,stuff)

stuff.ylabel='Sphericity';
stuff.save_name=[save_dir,'/Sph-T'];
stuff.color='blue';
plot_generic(numbers,sphericity,stuff)

stuff.ylabel='CC_{r}';
stuff.save_name=[save_dir,'/nxcc-T'];
stuff.color='red';
plot_generic(numbers,nxcc,stuff)

stuff.ylabel='CC_{d}';
stuff.save_name=[save_dir,'/dxcc-T'];
stuff.color='green';
if no_data == 0, plot_generic(numbers,dxcc,stuff),end

%stuff.ylabel='Amplitude Spread (\sigma / \mu)';
%stuff.save_name=[save_dir,'/amp-spread'];
%stuff.color='b';
%plot_generic(numbers,std_dev./amp_mean,stuff)

stuff.ylabel='Phase skewness';
stuff.save_name=[save_dir,'/ph-skew'];
stuff.color='b';
plot_generic(numbers,skew_ph,stuff)

%stuff.ylabel='Phase (Avg. per pixel)';
%stuff.save_name=[save_dir,'/strain-T'];
%stuff.color='r';
%plot_generic(numbers,strain,stuff)

stuff.ylabel='Phase (radians)';
stuff.save_name=[save_dir,'/strain-1sig-T'];
stuff.color='b';
plot_generic(numbers,strain_sig,stuff)

% stuff.ylabel='Phase (radians)';
% stuff.save_name=[save_dir,'/strain-2-T'];
% stuff.color='r';
% stuff.color2='b';
% stuff.legend={'\mu','\sigma'};
% plot_generic(numbers,strain,stuff,strain_sig)

stuff.ylabel='Corr. Coeff.';
stuff.save_name=[save_dir,'/CC-2-T'];
stuff.color='g';
stuff.color2='r';
stuff.legend={'CC_{d}','CC_{r}'};
if no_data == 0, plot_generic(numbers,dxcc,stuff,nxcc),end

stuff.ylabel='Corr. Coeff./\Delta V';
stuff.save_name=[save_dir,'/CC-2-V-T'];
stuff.color='g';
stuff.color2='r';
stuff.color3='b';
stuff.legend={'CC_{d}','CC_{r}','\Delta V'};
if no_data == 0,plot_generic(numbers,dxcc,stuff,nxcc,volume),end

if do_point == 1
    stuff.ylabel='Phase';
    stuff.save_name=[save_dir,'/ph-point'];
    stuff.color='b';
    plot_generic(numbers,angle(track_point),stuff)
end
%% plot the histograms as lines only
for qq=1:numel(amp_fs)
    
    hh_ph=atan2(imag(temp_series(:,:,:,qq)),real(temp_series(:,:,:,qq)));
    hh_qq=hh_ph(find(vol_series(:,:,:,qq) ~= 0) );
    [h xh]=imhist( (hh_qq(:)+pi)/2/pi,ph_bins);
    hist_ph(:,qq)=h(:);
    
    hh_amp=abs(temp_series(:,:,:,qq));
    hh_qq=hh_amp(find(vol_series(:,:,:,qq) ~= 0) );
    [h xl]=imhist( (hh_qq(:)),50);
    hist_amp(:,qq)=h(:);
    
    %fit a gauss to the histogram
    [gg x] = fit_gauss_data(hist_amp(10:end,qq));
    full=[zeros([1,9]),gg];
    full(full == 0)=NaN;
    hist_fit_amp(:,qq)=full;
    hist_fit_params(:,qq)=x;
end

stuff.legend=[];
stuff.ylabel='Amplitude Spread (\sigma / \mu)';
stuff.save_name=[save_dir,'/amp-spread'];
stuff.color='b';
plot_generic(numbers,hist_fit_params(2,:)./hist_fit_params(1,:),stuff)

if hist_range ~= 1
    nhist=ph_bins;
    %hist_amp_temp=hist_amp(nhist*(.5-hist_range/2):nhist*(.5+hist_range/2),:);

    hist_ph_temp=hist_ph(nhist*(.5-hist_range/2):nhist*(.5+hist_range/2),:);
    
    xh_temp=xh(nhist*(.5-hist_range/2):nhist*(.5+hist_range/2));
    %xl_temp=xl(nhist*(.5-hist_range/2):nhist*(.5+hist_range/2));
    
    xh=xh_temp;
    %xl=xl_temp;
    
    %hist_amp=hist_amp_temp;
    hist_ph=hist_ph_temp;
    
end

name=[save_dir,'/Ph-hist'];
legend_str=numbers_string;

if offset == 1,
    name=[name,'-OS'];
    hist_ph=add_offset_line_series(hist_ph,numel(amp_fs),offset_ph*max(hist_ph(:)),1);end

ph_max=max([max(hist_ph(:))]);
ph_min=min([min(hist_ph(:))]);
axis_r=[min(xh) max(xh) 1.1*ph_min 1.1*ph_max]; 
plot_amp_mult_color(hist_ph,1*size(hist_ph,2),lw,font_size,name,axis_r,legend_str,'Frequency',xh,'Phase (2 \pi)',stuff)

name=[save_dir,'/Amp-hist'];

if offset == 1,
    name=[name,'-OS'];
    hist_amp_old=hist_amp;
    hist_amp=add_offset_line_series(hist_amp,numel(amp_fs),offset_amp*max(hist_amp(:)),1);end

a_max=max([max(hist_amp(:))]);
a_min=min([min(hist_amp(:))]);
axis_r=[min(xl) max(xl) 1.1*a_min 1.1*a_max]; 

plot_amp_mult_color(hist_amp,1*size(hist_amp,2),lw,font_size,name,axis_r,[],'Frequency',xl,'Amplitude',stuff)

%plot amp hist with fitted curves
hist_amp=hist_amp_old;
name=[save_dir,'/Amp-hist-fit'];
if offset == 1,
    name=[name,'-OS'];
    hist_amp=add_offset_line_series(hist_amp,numel(amp_fs),offset_amp*max(hist_amp(:)),1);
    
    for qq=1:n_temps
        hist_fit_amp(:,qq)=hist_fit_amp(:,qq)-max(hist_fit_amp(:,qq))+max(hist_amp(:,qq));
    end
end

ww=1;
vv=1;
for qq=1:2*n_temps
    if mod(qq,2) == 1
        hist_amp_and_fit(:,qq)=hist_amp(:,ww);
        ww=ww+1;
    else
        hist_amp_and_fit(:,qq)=hist_fit_amp(:,vv);
        vv=vv+1;
    end
end
        
a_max=max([max(hist_amp_and_fit(:))]);
a_min=min([min(hist_amp_and_fit(:))]);
axis_r=[min(xl) max(xl) 1.1*a_min 1.1*a_max]; 

plot_amp_mult_color(hist_amp_and_fit,1*size(hist_amp_and_fit,2),lw,font_size,name,axis_r,[],'Frequency',xl,'Amplitude',stuff)

%


%% plot the line-outs

name=[save_dir,'/LO-Amp-X'];
if offset == 1,
    name=[name,'-OS'];
    x_amp_lines=add_offset_line_series(x_amp_lines,numel(amp_fs),.3,1);end

axis_r=[xls/2-xls/4 xls/2+xls/4 0 1.1*max(x_amp_lines(:))]; 
legend_str=numbers_string;
plot_amp_mult_color(x_amp_lines,1*size(x_amp_lines,2),lw,font_size,name,axis_r,legend_str,'Amplitude (Arb. units)',[],[],stuff)


name=[save_dir,'/LO-Amp-Y'];
if offset == 1,
    name=[name,'-OS'];
    y_amp_lines=add_offset_line_series(y_amp_lines,numel(amp_fs),.3,1);end
axis_r=[yls/2-yls/4 yls/2+yls/4 0 1.1*max(y_amp_lines(:))]; 
legend_str=numbers_string;
plot_amp_mult_color(y_amp_lines,1*size(y_amp_lines,2),lw,font_size,name,axis_r,[],'Amplitude (Arb. units)',[],[],stuff)


name=[save_dir,'/LO-Amp-Z'];
if offset == 1,
    name=[name,'-OS'];
    z_amp_lines=add_offset_line_series(z_amp_lines,numel(amp_fs),.3,1);end
axis_r=[zls/2-zls/4 zls/2+zls/4 0 1.1*max(z_amp_lines(:))]; 
legend_str=numbers_string;
plot_amp_mult_color(z_amp_lines,1*size(z_amp_lines,2),lw,font_size,name,axis_r,[],'Amplitude (Arb. units)',[],[],stuff)

% do the phases
ph_max=max([max(x_ph_lines(:)),max(y_ph_lines(:)),max(z_ph_lines(:))]);
ph_min=min([min(x_ph_lines(:)),min(y_ph_lines(:)),min(z_ph_lines(:))]);

name=[save_dir,'/LO-Ph-X'];
if offset == 1,
    name=[name,'-OS'];
    x_ph_lines=add_offset_line_series(x_ph_lines,numel(amp_fs),.5,1);
    ph_max=max(x_ph_lines(:));ph_min=min(x_ph_lines(:));
end
axis_r=[xls/2-xls/4 xls/2+xls/4 1.1*ph_min 1.1*ph_max]; 
legend_str=numbers_string;
plot_amp_mult_color(x_ph_lines,1*size(x_ph_lines,2),lw,font_size,name,axis_r,legend_str,'Phase (radians)',[],[],stuff)

name=[save_dir,'/LO-Ph-Y'];
if offset == 1,
    name=[name,'-OS'];
    y_ph_lines=add_offset_line_series(y_ph_lines,numel(amp_fs),.5,1);
    ph_max=max(y_ph_lines(:));ph_min=min(y_ph_lines(:));
end
axis_r=[yls/2-yls/4 yls/2+yls/4 1.1*ph_min 1.1*ph_max]; 
legend_str=numbers_string;
plot_amp_mult_color(y_ph_lines,1*size(y_ph_lines,2),lw,font_size,name,axis_r,[],'Phase (radians)',[],[],stuff)

name=[save_dir,'/LO-Ph-Z'];
if offset == 1,
    name=[name,'-OS'];
    z_ph_lines=add_offset_line_series(z_ph_lines,numel(amp_fs),.6,1);
    ph_max=max(z_ph_lines(:));ph_min=min(z_ph_lines(:));
end
axis_r=[zls/2-zls/4 zls/2+zls/4 1.1*ph_min 1.1*ph_max]; 
legend_str=numbers_string;
plot_amp_mult_color(z_ph_lines,1*size(z_ph_lines,2),lw,font_size,name,axis_r,[],'Phase (radians)',[],[],stuff)


% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% plot(x_amp_lines,'LineWidth',lw)
% set(gca,'FontSize',round(0.8*font_size))
% xlabel('Pixel','FontSize', font_size), ylabel('Amplitude (Arb. units)','FontSize', font_size)
% axis([xls/2-xls/4 xls/2+xls/4 0 1.1 ])
% legend(numbers_string)
% saveas(fh, [save_dir,'/LO-Amp-X'], 'epsc');
% print(fh, '-dpng','-r600', [save_dir,'/LO-Amp-X']);
% 
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% plot(y_amp_lines,'LineWidth',lw)
% set(gca,'FontSize',round(0.8*font_size))
% xlabel('Pixel','FontSize', font_size), ylabel('Amplitude (Arb. units)','FontSize', font_size)
% axis([yls/2-yls/4 yls/2+yls/4 0 1.1 ])
% legend(numbers_string)
% saveas(fh, [save_dir,'/LO-Amp-Y'], 'epsc');
% print(fh, '-dpng','-r600', [save_dir,'/LO-Amp-Y']);
% 
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% plot(z_amp_lines,'LineWidth',lw)
% set(gca,'FontSize',round(0.8*font_size))
% xlabel('Pixel','FontSize', font_size), ylabel('Amplitude (Arb. units)','FontSize', font_size)
% axis([zls/2-zls/4 zls/2+zls/4 0 1.1 ])
% legend(numbers_string)
% saveas(fh, [save_dir,'/LO-Amp-Z'], 'epsc');
% print(fh, '-dpng','-r600', [save_dir,'/LO-Amp-Z']);



% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% plot(x_ph_lines,'LineWidth',lw)
% set(gca,'FontSize',round(0.8*font_size))
% xlabel('Pixel','FontSize', font_size), ylabel('Phase (radians)','FontSize', font_size)
% axis([xls/2-xls/4 xls/2+xls/4 1.1*ph_min 1.1*ph_max ])
% legend(numbers_string)
% saveas(fh, [save_dir,'/LO-Ph-X'], 'epsc');
% print(fh, '-dpng','-r600', [save_dir,'/LO-Ph-X']);
% 
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% plot(y_ph_lines,'LineWidth',lw)
% set(gca,'FontSize',round(0.8*font_size))
% xlabel('Pixel','FontSize', font_size), ylabel('Phase (radians)','FontSize', font_size)
% axis([yls/2-yls/4 yls/2+yls/4 1.1*ph_min 1.1*ph_max ])
% legend(numbers_string)
% saveas(fh, [save_dir,'/LO-Ph-Y'], 'epsc');
% print(fh, '-dpng','-r600', [save_dir,'/LO-Ph-Y']);
% 
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% plot(z_ph_lines,'LineWidth',lw)
% set(gca,'FontSize',round(0.8*font_size))
% xlabel('Pixel','FontSize', font_size), ylabel('Phase (radians)','FontSize', font_size)
% axis([zls/2-zls/4 zls/2+zls/4 1.1*ph_min 1.1*ph_max ])
% legend(numbers_string)
% saveas(fh, [save_dir,'/LO-Ph-Z'], 'epsc');
% print(fh, '-dpng','-r600', [save_dir,'/LO-Ph-Z']);
close all

end


function pn = load_rec(amp_f,ph_f)

load(amp_f,'-mat')
amp=array;
array=0;
load(ph_f,'-mat');
ph=array;
array=0;

pn=amp.*exp(i*ph);

end

function [shifted]=sub_pixel_shift(array,row_shift,col_shift,z_shift)

buf2ft=fftn(array);
[nr,nc,nz]=size(buf2ft);
Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
Nz = ifftshift([-fix(nz/2):ceil(nz/2)-1]);
[Nc,Nr,Nz] = meshgrid(Nc,Nr,Nz);
Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc-z_shift*Nz/nz));
%Greg = Greg*exp(i*diffphase);
shifted=ifftn(Greg);

end

function cnj = conj_reflect(array)

F=ifftshift(fftn(fftshift(array)));

cnj=ifftshift(ifftn(fftshift(conj(F))));


end

function cnj_rf=is_conj_ref(a,b)

c1=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(conj_reflect(b),.1,.1));
c2=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(b,.1,.1));

c1=max(c1(:));
c2=max(c2(:));

if c1 > c2,cnj_rf=1;else cnj_rf=0;end


end

function plot_generic(x,y,stuff,y2,y3)

AA=blue2red;
AA=circshift(AA,[64,0]);
BB=reverse(AA(16+1:128-16,:));
ind=round(1:(size(BB,1)/(8)):size(BB,1));
ColorSet=BB(ind,:);

try
    y2;
catch
    y2=[];
end
try
    y3;
catch
    y3=[];
end

try
    stuff.legend;
catch
    stuff.legend=[];
end

try
    lw=stuff.lw;
catch
    lw=1.5;
end
try
    font_size=stuff.font_size;
catch
    font_size=25;
end
try
    xlab=stuff.xlabel;
catch
    xlab='';
end
try
    ylab=stuff.ylabel;
catch
    ylab='';
end
try
    save_name=stuff.save_name;
catch
    save_name=[];
end

try
    stuff.color;
catch
    stuff.color='blue';
end
try
    stuff.color2;
catch
    stuff.color2='red';
end
try
    stuff.color3;
catch
    stuff.color3='green';
end
try
    stuff.n_temps_cold;
catch
    stuff.n_temps_cold=0;
end

try
    stuff.dots_only;  %plot data points only
catch
    stuff.dots_only=0;
end
try
    stuff.bold_plot;    %plot in bold
catch
    stuff.bold_plot=1;
end
%
markersize=10;
%
fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

if stuff.n_temps_cold == 0
    if stuff.dots_only == 0
        plot(x,y,'LineWidth',lw,'Color',stuff.color)   
    else
        plot(x,y,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color,'MarkerEdgeColor',stuff.color)
    end
else
    if stuff.dots_only == 0
        plot(x(1:end-stuff.n_temps_cold),y(1:end-stuff.n_temps_cold),'LineWidth',lw,'Color',ColorSet(4,:))
        hold on
        plot(x(end-stuff.n_temps_cold:end),y(end-stuff.n_temps_cold:end),'LineWidth',lw,'Color',ColorSet(1,:))
    else
        plot(x(1:end-stuff.n_temps_cold),y(1:end-stuff.n_temps_cold),'o','MarkerSize',8,'MarkerFaceColor',ColorSet(4,:))
        hold on
        plot(x(end-stuff.n_temps_cold:end),y(end-stuff.n_temps_cold:end),'o','MarkerSize',8,'MarkerFaceColor',ColorSet(1,:))
    end
end

set(gca,'FontSize',round(0.8*font_size))
if stuff.bold_plot == 1
    xlabel(xlab,'FontSize', font_size,'FontWeight','bold');ylabel(ylab,'FontSize', font_size,'FontWeight','bold');
else
    xlabel(xlab,'FontSize', font_size);ylabel(ylab,'FontSize', font_size); 
end

if numel(y2) ~= 0
    hold on;
    if stuff.dots_only == 0
        plot(x,y2,'LineWidth',lw,'Color',stuff.color2)
    else
        plot(x,y2,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color2,'MarkerEdgeColor',stuff.color2)
    end
end

if numel(y3) ~= 0
    if stuff.dots_only == 0
        plot(x,y3,'LineWidth',lw,'Color',stuff.color3)
    else
        plot(x,y3,'o','MarkerSize',markersize,'MarkerFaceColor',stuff.color3,'MarkerEdgeColor',stuff.color3)
    end
end
hold off;

[aa]=axis;
if numel(y2) == 0,y2=y;end
if numel(y3) == 0,y3=y;end

mm=min([min(y3(:)),min(y2(:)),min(y(:))]);
nn=max([max(y3(:)),max(y2(:)),max(y(:))]);

aa(3)=mm-.1*abs(mm-nn);
aa(4)=nn+.1*abs(mm-nn);

axis([aa(1) aa(2) aa(3) aa(4)])

if numel(stuff.legend) ~= 0
    legend(stuff.legend,'Location','best')
end
box off

if stuff.bold_plot == 1,set(gca,'FontWeight','bold','LineWidth',.8*lw);end
%else 
    %set(gca,'LineWidth',.8*lw);end

if numel(save_name) ~= 0
    saveas(fh, save_name, 'epsc');
    print(fh, '-dpng','-r600', save_name);
end



end

function plot_amp_mult_color(data,n_colors,lw,font_size,name,axis_r,legend_str,ylab,x,xlab,stuff)

try
    x;
catch
    x=[];
end
try
    xlab;
    if isempty(xlab), xlab='Pixel';end
catch
    xlab='Pixel';
end


try
    offset;
catch
    offset=0;
end
try
    stuff.bold_plot;    %plot in bold
catch
    stuff.bold_plot=1;
end

n_plots=size(data,2);

AA=blue2red;
AA=circshift(AA,[64,0]);
BB=reverse(AA(16+1:128-16,:));
ind=round(1:(size(BB,1)/(n_plots)):size(BB,1));
ColorSet=BB(ind,:);
%ColorSet=(varycolor(n_colors));
%AA=reverse(ColorSet(1:n_plots,:));
%ColorSet=0;
%ColorSet=AA;

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 

hold on

if numel(x) ~= 0
    plot(x,data(:,1),'LineWidth',lw,'Color',ColorSet(1,:))
else
    plot(data(:,1),'LineWidth',lw,'Color',ColorSet(1,:))
end

if stuff.bold_plot==1
    set(gca,'FontSize',round(0.8*font_size),'FontWeight','bold','LineWidth',.8*lw)
    xlabel(xlab,'FontSize', font_size,'FontWeight','bold'), ylabel(ylab,'FontSize', font_size,'FontWeight','bold')
else
    set(gca,'FontSize',round(0.8*font_size))
    xlabel(xlab,'FontSize', font_size), ylabel(ylab,'FontSize', font_size)        
end

axis(axis_r)

for qq=2:n_plots
    if numel(x) ~= 0
        plot(x,data(:,qq),'LineWidth',lw,'Color',ColorSet(qq,:))
    else
        plot(data(:,qq),'LineWidth',lw,'Color',ColorSet(qq,:))
    end
end
hold off

%set(gca,'FontSize',round(0.8*font_size))
%xlabel('Pixel','FontSize', font_size), ylabel('Amplitude (Arb. units)','FontSize', font_size)
%axis(axis_r)
if numel(legend_str) ~= 0
    legend(legend_str,'Location','BestOutside')
end
box off
%set(gca,'FontWeight','bold','LineWidth',.8*lw)


saveas(fh, name, 'epsc');
print(fh, '-dpng','-r600', name);

end

function out=add_offset_line_series(data,plots,offset,reverse)

try
    reverse;
catch
    reverse=0;
end

if reverse == 0,
    for qq=1:plots
        out(:,qq)=data(:,qq)+(qq-1)*offset;
    end
else
    for qq=1:plots
        out(:,qq)=data(:,qq)+(plots-qq)*offset;
    end
end

end

function output_file_names(save_name,files)

n_files=max(size(files));

fid=fopen([save_name],'w');

for kk=1:n_files

    fprintf(fid,'\n');
    fprintf(fid,char(files(kk).name));
    fprintf(fid,'\n');
    
end

fclose(fid)

end

function output_string(save_name,string)

fid=fopen([save_name],'w');

fprintf(fid,'\n');
fprintf(fid,'%c',string);
fprintf(fid,'\n');

fclose(fid)

end

function sz = plot_amp_phase(pn,save_dir,ph_range,sz)

if ischar(save_dir) == 0
   
    fdir=strtrim(char(save_dir(1)));
    sdir=strtrim(char(save_dir(2)));
    tdir=strtrim(char(save_dir(3)));
    
else
    
    fdir=save_dir;
    sdir=save_dir;
    tdir=save_dir;
    
end

%

nx=size(pn,2);
ny=size(pn,1);
nz=size(pn,3);

switch sum(sz(:))
    
    case 0
        fact=1.5;
        num_x=abs(pn(round(ny/2),:,round(nz/2)));
        num_y=abs(pn(:,round(nx/2),round(nz/2)));
        num_z=abs(pn(round(ny/2),round(nx/2),:));
        ind_x=(num_x >= .1*max(num_x));
        ind_y=(num_y >= .1*max(num_y));
        ind_z=(num_z >= .1*max(num_z));
        sx=fact*sum(ind_x(:));
        sy=fact*sum(ind_y(:));
        sz=fact*sum(ind_z(:));

        xx=round([nx/2-sx/2,nx/2+sx/2]);
        xx(xx < 1)=1;
        xx(xx > nx)=nx;
        yy=round([ny/2-sy/2,ny/2+sy/2]);
        yy(yy < 1)=1;
        yy(yy > ny)=ny;
        zz=round([nz/2-sz/2,nz/2+sz/2]);
        zz(zz < 1)=1;
        zz(zz > nz)=nz;
        sz=[xx,yy,zz];
    case 1
        xx=[1,nx];%[nx/2-nx/3,nx/2+nx/3];
        yy=[1,ny];%[ny/2-ny/3,ny/2+ny/3];
        zz=[1,nz];%[nz/2-nz/3,nz/2+nz/3];
    otherwise
        xx=[sz(1),sz(2)];
        yy=[sz(3),sz(4)];
        zz=[sz(5),sz(6)];
end
    

lw=1.5;

phase=atan2(imag(pn),real(pn) );
amp=abs(pn);


fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
%ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
ph = extract_3D_slice(amp,'xy' );
imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
axis equal
%caxis(ph_range)
saveas(fh, [fdir,'Amp-xy'],'epsc');
print(fh, '-dpng','-r300', [fdir,'Amp-xy']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(amp,'xz' );
imagesc(ph(zz(1):zz(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
axis equal
%caxis(ph_range)
saveas(fh, [sdir,'Amp-xz'],'epsc');
print(fh, '-dpng','-r300', [sdir,'Amp-xz']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(amp,'zy' );
imagesc(ph(yy(1):yy(2),zz(1):zz(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
axis equal
%caxis(ph_range)
saveas(fh, [tdir,'Amp-zy'],'epsc');
print(fh, '-dpng','-r300', [tdir,'Amp-zy']);

if isreal(pn) ~= 1
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    %ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
    ph = extract_3D_slice(phase,'xy' );
    imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
    h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
    axis equal
    caxis(ph_range);
    saveas(fh, [fdir,'Ph-xy'],'epsc');
    print(fh, '-dpng','-r300', [fdir,'Ph-xy']);

    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    ph = extract_3D_slice(phase,'xz' );
    imagesc(ph(zz(1):zz(2),xx(1):xx(2)));
    h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
    axis equal
    caxis(ph_range);
    saveas(fh, [sdir,'Ph-xz'],'epsc');
    print(fh, '-dpng','-r300', [sdir,'Ph-xz']);

    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    ph = extract_3D_slice(phase,'zy' );
    imagesc(ph(yy(1):yy(2),zz(1):zz(2)));
    h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
    axis equal
    caxis(ph_range);
    saveas(fh, [tdir,'Ph-zy'],'epsc');
    print(fh, '-dpng','-r300', [tdir,'Ph-zy']);
end
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% imagesc(abs(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))))
% h=colorbar('location','EastOutside','fontsize',0,'fontweight','bold');
% saveas(fh, [save_dir,'-Amp-c'],'epsc');
% print(fh, '-dpng','-r300', [save_dir,'-Amp-c']);


end
function [gg x] = fit_gauss_data(data)

options = optimset('Display','off','Algorithm','interior-point');
 
tempx=find(data == max(data));
x0(1)=tempx(1);

x0(2)=2;

x0(3)=max(data);

x0(4)=data(1);
%x0(5)=1/20*x0(3);

lb=[0,0,0,-x0(3)];%,-5,-x0(3)];
ub=[numel(data),numel(data),5*x0(3),2*x0(3)];%,5,x0(3)];

f=@(x)gauss_fit(x,data);
x=fmincon(f,x0,[],[],[],[],[lb],[ub],[],options);


mean=x(1);
sig=x(2);
A0=x(3);
%m=x(4);
c=x(4);

xx=1:numel(data);   %abisca values


gg=A0*exp(-0.5*(xx-mean).^2/sig.^2)+c;


end
function E = gauss_fit(x,data,range )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%y=A0*exp(-0.5 x*x/(sig*sig))+mx+c

mean=x(1);
sig=x(2);
A0=x(3);
%m=x(4);
c=x(4);


xx=1:numel(data);   %abisca values


gauss=A0*exp(-0.5*(xx-mean).^2/sig.^2)+c;

data=data(:);
gauss=gauss(:);

E=sum(abs(gauss-data));


end
% function plot_amp_phase(pn,save_dir,ph_range)
% 
% nx=size(pn,2);
% ny=size(pn,1);
% nz=size(pn,3);
% 
% xx=[1,nx];%[.15,.85]*nx;%[nx/2-nx/4,nx/2+nx/4];
% yy=[1,ny];%[.15,.85]*ny;%[ny/2-ny/4,ny/2+ny/4];
% zz=[1,nz];%[.15,.85]*nz;%[nz/2-nz/4,nz/2+nz/4];
% 
% lw=1.5;
% 
% phase=atan2(imag(pn),real(pn) );
% 
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% ph=atan2(imag(pn(yy(1):yy(2),xx(1):xx(2),round(nz/2)),real(pn(yy(1):yy(2),xx(1):xx(2),round(nz/2))) );
% 
% imagesc(ph)
% h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
% caxis(ph_range)
% saveas(fh, [save_dir,'-Ph-xy-c'],'epsc');
% print(fh, '-dpng','-r300', [save_dir,'-Ph-xy-c']);
% 
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% ph = extract_3D_slice(phase,'xz' );
% imagesc(ph(zz(1):zz(2),xx(1):xx(2)));
% h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
% caxis(ph_range)
% saveas(fh, [save_dir,'-Ph-xz-c'],'epsc');
% print(fh, '-dpng','-r300', [save_dir,'-Ph-xz-c']);
% 
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% ph = extract_3D_slice(phase,'yz' );
% imagesc(ph(zz(1):zz(2),yy(1):yy(2)));
% h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
% caxis(ph_range)
% saveas(fh, [save_dir,'-Ph-yz-c'],'epsc');
% print(fh, '-dpng','-r300', [save_dir,'-Ph-yz-c']);
% 
% 
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% imagesc(abs(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))))
% h=colorbar('location','EastOutside','fontsize',0,'fontweight','bold');
% saveas(fh, [save_dir,'-Amp-c'],'epsc');
% print(fh, '-dpng','-r300', [save_dir,'-Amp-c']);
% 
% 
% end