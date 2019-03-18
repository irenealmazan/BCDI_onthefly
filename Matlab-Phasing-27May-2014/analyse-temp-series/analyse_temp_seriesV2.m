function [ output_args ] = analyse_temp_seriesV2(top_dir,prefix,params)

%jclark
% to analyse a temperature series of CDI
% it is assumed the reconstruction to be analysed will be an average
% from several random starts.  with the summation performed using
% sum_align_reconstructions
% top_dir should be the dir that has the temperatures below it
% e.g. top_dir/T030/rand-starts/CVL/

params=set_params_defaults(params);

no_data=params.no_data;    %don't do the data, its takes time, use data = 0
avgr=params.avgr;       %region to average over, nxn
offset=params.offset;     %plot the lineouts with an offset?, yes=1
shell_hist=params.shell_hist;   %do the shell histograms, yes=1
amp_norm=params.amp_norm;       %normalise amplitudes of reconstruction to that of the first
                                %such that max(first)=1, and then
                                %sum(second)=sum(first).....
plot_cold=params.plot_cold;    %set =1 to plot the cooled down temperature 30 if present
cold_string=params.cold_string;  %s string that is in the file name of the cold temps that distniguishes them from the rest
lab=params.lab;          %load the transformed reconstruction
vtk_out=params.vtk_out;      %output vtk file (will be correctly orientated)
auto_size_plot=params.auto_size_plot;   %auto size (=0), use all =1
conj_first=params.conj_first;     %conjugate the reference? (=1 yes, 0=no)
rem_ramp=params.rem_ramp;     %remove ramp in sample plane, =1 yes
upsize=params.upsize;       %upsample factor
z_phase=params.z_phase;      %set mean to zero (1=yes)
apply_sup=params.apply_sup;    %multiply by support (vtk output)
hist_range=params.hist_range;     %fraction of 2 pi
ph_bins=params.ph_bins;
ph_range=params.ph_range;  %range for output image
offset_amp=params.offset_amp;          %offset to plot amp histogam as fraction of max
offset_ph=params.offset_ph;          %offset to plot ph histogam as fraction of max
vsig=params.vsig;               %calc values for the volume
vth=params.vth;
stuff.bold_plot=params.stuff.bold_plot;
stuff.dots_only=params.stuff.dots_only;
sort_numerically=params.sort_numerically;     %will sort according to number.  do not use if heating series has coold down temps, will lose the order
nnumbs_out=params.nnumbs_out;%[2];       %number of digits in the numbers for save name output, use [] to turn off
string_prefix=params.string_prefix;
custom_xaxis=params.custom_xaxis;        %use a custom axis, leave as [] for nothing
cxlabel=params.cxlabel;%'
points=params.points;
looffset=params.line_out_cent_offset;
outdir=params.outdir;
calc_diff=params.calc_diff;
apply_sup_slices=params.apply_sup_slices;
interp_slices=params.interp_slices;
upsamp_vtk=params.upsamp_vtk;
box_vtk=params.box_vtk;
nshells=params.nshells;                %number of shells to use
params.maxy=params.shell_maxy;


params_in=params;

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
    params.maxy;
catch
    params.maxy=1000;
end
maxy=params.maxy;

try
    top_dir;
catch
    name_of_this_file='analyse_temp_seriesV2';
    dir_file=which(name_of_this_file);    %finds the current folder for the phasing
    top_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
end
try
    nnumbs_out;
catch
    nnumbs_out=[];
end

try
    params.cold_string;
    cold_string=params.cold_string;
catch
    cold_string='a/R';
end
try
    plot_subset=params.plot_subset;
catch
    plot_subset=[];
end

if auto_size_plot == 1,box_vtk=0;end

% sort the names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% get names of the summed and averaged reconstructions and a paramter file

if lab ~= 1,string_amp=[prefix,'**-AMP.rec'];else string_amp=[prefix,'**LAB.rec'];end
amp_fs=rdir([top_dir,'**/*/',string_amp]);  %amp file names
n_temps=size(amp_fs,1);                         %get the number of temperatures

string_ph=[prefix,'**-PH.rec'];
ph_fs=rdir([top_dir,'**/*/',string_ph]);    %ph file names

amp_fs_c=rdir([top_dir,string_prefix,'*',cold_string,'*/*/',string_amp]); %get the names of the cold ones, if any
n_temps_cold=size(amp_fs_c,1);                 %get the number of cold temps, an use this as a sanity check later on
check=zeros([1,n_temps]);

string_params=[prefix,'**PARAMS.mat'];      %params file names
pms_fs=rdir([top_dir,'**/*/',string_params]);

%check which names match the 'cooling down' string, these will go at the
%end
for qq=1:n_temps
    if numel(regexp(char(amp_fs(qq).name),cold_string,'match')) > 0,check(qq)=1;end
    %careful of the lower or upper case a/r or a/R
end

%check if the number of cold names matches from two different methods
if sum(check) ~= n_temps_cold,disp('ERROR determining cold temperatures....'), else disp('Cold temperatures determined....'),end

%%

warm_ind=find(check == 0);   %the indexes of the not cold ones
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

% start loading reconstructions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
surface_area=0;
sphericity=0;
nxcc=0;
strain=0;
strain_sig=0;
std_dev=0;
numbers=0;
skew_ph=0;
numbers_string={};
params.colors=['r','b'];        %used for the histogram colors
%nshells=2;                      %number of shells for the histogram

save_dir=[top_dir,outdir];
vtk_dir=[save_dir,'/vtks'];
make_directory({save_dir,vtk_dir});

if calc_diff == 1,
    vtk_diff_dir=[save_dir,'/diff/vtks'];
    diff_dir_xy=[save_dir,'/diff/slices/xy'];diff_dir_xz=[save_dir,'/diff/slices/xz'];diff_dir_zy=[save_dir,'/diff/slices/zy'];          
    make_directory({vtk_diff_dir,diff_dir_xy,diff_dir_xz,diff_dir_zy});   
end

xy_dir=[save_dir,'/xy'];xz_dir=[save_dir,'/xz'];zy_dir=[save_dir,'/zy'];
%dir for the linout locations
LO_dir=[save_dir,'/LO'];
make_directory({xy_dir,xz_dir,zy_dir,LO_dir});


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
surface_area(1)=sum(sum(sum(bwperim(vol_sup))));

%output a vtk into the dir if specified
if vtk_out == 1
    vtkname=[vtk_dir,'/',string_prefix,nums,'-Amp-Phase.vtk'];   
    if apply_sup == 1,vtkarray=flip_all_dim(first.*vol_sup);disp('Applying support to vtk output....')
    else vtkarray=flip_all_dim(first);end
    if upsamp_vtk ~= 1,vtkarray=FFT_interpolate(vtkarray,round(upsamp_vtk*[size(vtkarray,2),size(vtkarray,1),size(vtkarray,3)]));vtkarray=vtkarray.*shrink_wrap(abs(vtkarray),.05,.5);end   
end

if apply_sup_slices == 1,first=first.*vol_sup;disp('Applying support to slices....');end %first(abs(first) == 0)=NaN;

%calc corr coeff for data
if no_data == 0, dxcc(1)=max(max(max(normxcorr3(data, data, 'same'))));end %normalised cross-corr data

% calc the 'strain'.  it is the integrated phase, would need to be
% unwrapped for large phase
ph=atan2(imag(first),real(first));
ph=ph-mean( ph( find(ph.*vol_sup ~=0)) );
strain(1)=sqrt(sum(sum(sum(abs(vol_sup.*ph).^2)))/numel(find(vol_sup > 0)));
strain_sig(1)=std( ph(find(ph.*vol_sup ~=0)));
skew_ph(1)=skewness( ph(find(ph.*vol_sup ~=0)));


%output phase and amplitude and volume support
save_dirs={[xy_dir,'/',nums,'I-'],[xz_dir,'/',nums,'I-'],[zy_dir,'/',nums,'I-']};
save_dirsS={[xy_dir,'/',nums,'S-'],[xz_dir,'/',nums,'S-'],[zy_dir,'/',nums,'S-']};
asp = plot_amp_phase(first,[save_dirs],ph_range,auto_size_plot,params_in);%.*vol_sup
plot_amp_phase(vol_sup,[save_dirsS],ph_range,asp,params_in);

%output vtk, needed to get asp before saving though
if box_vtk == 1,
    asp_a=asp*upsamp_vtk;
    box_vtk_s=zeros(size(vtkarray));box_vtk_s(asp_a(3):asp_a(4),asp_a(1):asp_a(2),asp_a(5):asp_a(6))=1;else box_vtk_s=1;end
if vtk_out == 1
    savevtk2scalar(box_vtk_s.*abs(vtkarray)/max(abs(vtkarray(:))),vtkname,angle(vtkarray),1)
    vtkarray=[];end

%output slices for a movie through the crystal
if params_in.output_ind_slices == 1
    nano_sl_dir=[save_dir,'/nanoslices/'];
    output_N_slices(angle(first),params_in.output_ind_N,'xy',asp,[nano_sl_dir,nums,'/xy/'],'Ph-xy',ph_range);
    output_N_slices(angle(first),params_in.output_ind_N,'xz',asp,[nano_sl_dir,nums,'/xz/'],'Ph-xz',ph_range);
    output_N_slices(angle(first),params_in.output_ind_N,'zy',asp,[nano_sl_dir,nums,'/zy/'],'Ph-zy',ph_range);
        
    %make_avi_from_images([nano_sl_dir,nums,'/xy/Ph*'],[nano_sl_dir,nums,'/xy/Ph-xy'],5,500)
end

%asp = plot_amp_phase(first,[save_dirs,'/',nums],ph_range,auto_size_plot);%.*vol_sup
%plot_amp_phase(vol_sup,[save_dirs,'/',nums,'S'],ph_range,asp);

%calc histogram for amp + phase
params_sh=params_in;
params_sh.prefix='-A-';
%params_sh.colors=['r','b'];
params_sh.xaxis='Amplitude Value';
if shell_hist==1;
    params_sh=density_shells_histogram(abs(first),vol_sup,nshells,[save_dir,'/',string_prefix,char(numbers_string(1)),'/'],params_sh);
    %keep shell vals
    hparams.amp_shell_vals{1}=params_sh.shell_vals;
    params_sh.shell_vals=[];
end

params_sh.prefix='-P-';
%params_sh.colors=['c','b'];        %used for the histogram colors
params_sh.xaxis='Phase (2 \pi)';
if shell_hist == 1;
    params_sh=density_shells_histogram((ph+pi)/2/pi,vol_sup,nshells,[save_dir,'/',string_prefix,char(numbers_string(1)),'/'],params_sh);
    %keep shell vals
    hparams.ph_shell_vals{1}=params_sh.shell_vals;
    params_sh.shell_vals=[];
end

%calc some line-outs
if isempty(looffset) == 1,looffset=[0,0,0];end
yl=round(size(first,1)/2-(avgr-1)/2+looffset(2));
yh=round(size(first,1)/2+(avgr-1)/2+looffset(2));
zl=round(size(first,3)/2-(avgr-1)/2+looffset(3));
zh=round(size(first,3)/2+(avgr-1)/2+looffset(3));
xl=round(size(first,2)/2-(avgr-1)/2+looffset(1));
xh=round(size(first,2)/2+(avgr-1)/2+looffset(1));
    
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

%output a pic to show where the lineout is
save_dirsLO={[LO_dir,'/',nums,'LO-'],[LO_dir,'/',nums,'LO-'],[LO_dir,'/',nums,'LO-']};
amp_temp=abs(first);ph_temp=angle(first);
amp_temp(yl:yh,:,zl:zh)=max(amp_temp(:));amp_temp(:,xl:xh,zl:zh)=max(amp_temp(:));amp_temp(yl:yh,xl:xh,:)=max(amp_temp(:));
ph_temp(yl:yh,:,zl:zh)=max(ph_temp(:));ph_temp(:,xl:xh,zl:zh)=max(ph_temp(:));ph_temp(yl:yh,xl:xh,:)=max(ph_temp(:));
disp('Saving image with line-out locations....')

if params_in.interp_slices == 1
    plot_amp_phase(amp_temp.*exp(i*ph_temp),[save_dirsLO],ph_range,asp,params_in);amp_temp=[];ph_temp=[];LO_arr=zeros(size(first));
else
    plot_amp_phase(amp_temp.*exp(i*ph_temp),[save_dirsLO],ph_range,asp,params_in);amp_temp=[];ph_temp=[];LO_arr=zeros(size(first));
end

LO_arr(yl:yh,:,zl:zh)=1;LO_arr(:,xl:xh,zl:zh)=1;LO_arr(yl:yh,xl:xh,:)=1;
LO_arr=flip_all_dim(LO_arr);savevtk2scalar(LO_arr,[LO_dir,'/',nums,'LO.vtk'])

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
    surface_area(qq)=sum(sum(sum(bwperim(vol_sup))));
    
    %calc some phase stuff, set mean to pi
    ph=atan2(imag(crystal),real(crystal));
    ph=ph-mean( ph(find(ph.*vol_sup ~=0)));
    strain(qq)=sqrt(sum(sum(sum(abs(vol_sup.*ph).^2)))/numel(find(vol_sup > 0)));
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
        
        if apply_sup == 1,vtkarray=flip_all_dim(crystal.*vol_sup);disp('Applying support to vtk output....')
        else vtkarray=flip_all_dim(crystal);end
        if upsamp_vtk ~= 1,vtkarray=FFT_interpolate(vtkarray,round(upsamp_vtk*[size(vtkarray,2),size(vtkarray,1),size(vtkarray,3)]));vtkarray=vtkarray.*shrink_wrap(abs(vtkarray),.05,.5);end
        
        savevtk2scalar(box_vtk_s.*abs(vtkarray)/max(abs(vtkarray(:))),vtkname,angle(vtkarray),1)
        
        vtkarray=[];
        
        if calc_diff == 1
            disp('Outputting difference between first and current....')
            vtkname=[vtk_diff_dir,'/',string_prefix,nums_a,'-D-Amp-Phase','.vtk'];
            %diff_c=(first-crystal).*vol_sup;
            %diff_c=(exp(i*angle(first))-exp(i*angle(crystal))).*vol_sup;
            
            %diff_c=exp(i*(angle(crystal)-angle(first))).*abs(crystal).*vol_sup;
            
            diff_c=exp(i*(angle(crystal)-angle(first))).*abs(first).*vol_series(:,:,:,1);
            
            
            
            vtkarray=flip_all_dim(zero_phase(diff_c));
            
            if isempty(params_in.crop_vtk) ~= 1
                vtkarray=zero_pad_ver3(vtkarray,params_in.crop_vtk(1),params_in.crop_vtk(2),params_in.crop_vtk(3));
            end
            
            savevtk2scalar(abs(vtkarray),vtkname,angle(vtkarray),1)
            vtkarray=[];
        end       
    end

    if apply_sup_slices == 1,crystal=crystal.*vol_sup;disp('Applying support to slices....');end %crystal(abs(crystal) == 0)=NaN;

    %output slices for a movie through the crystal
     if params_in.output_ind_slices == 1
        output_N_slices(angle(crystal),params_in.output_ind_N,'xy',asp,[nano_sl_dir,nums_a,'/xy/'],'Ph-xy',ph_range);
        output_N_slices(angle(crystal),params_in.output_ind_N,'xz',asp,[nano_sl_dir,nums_a,'/xz/'],'Ph-xz',ph_range);
        output_N_slices(angle(crystal),params_in.output_ind_N,'zy',asp,[nano_sl_dir,nums_a,'/zy/'],'Ph-zy',ph_range);
        %make_avi_from_images([nano_sl_dir,nums,'/xy/Ph*'],[nano_sl_dir,nums,'/xy/Ph-xy'],5,500)
    end
    
    %output the slcies
    save_dirs={[xy_dir,'/',nums_a,'I-'],[xz_dir,'/',nums_a,'I-'],[zy_dir,'/',nums_a,'I-']};
    save_dirsS={[xy_dir,'/',nums_a,'S-'],[xz_dir,'/',nums_a,'S-'],[zy_dir,'/',nums_a,'S-']};
    plot_amp_phase(crystal,[save_dirs],ph_range,asp,params_in);%.*vol_sup
    plot_amp_phase(vol_sup,[save_dirsS],ph_range,asp,params_in);
    
    %output slices of the difference if calc_diff == 1
    if exist('diff_c')
        slices_diff_dirs={[diff_dir_xy,'/',nums_a,'I-'],[diff_dir_xz,'/',nums_a,'I-'],[diff_dir_zy,'/',nums_a,'I-']};
        
        if params_in.calc_diff_slices == 1
            output_N_slices(angle(diff_c),params_in.output_ind_N,'xy',asp,[diff_dir_xy,nums_a,'/xy/'],'Ph-xy',ph_range);
            output_N_slices(angle(diff_c),params_in.output_ind_N,'xz',asp,[diff_dir_xz,nums_a,'/xz/'],'Ph-xz',ph_range);
            output_N_slices(angle(diff_c),params_in.output_ind_N,'zy',asp,[diff_dir_zy,nums_a,'/zy/'],'Ph-zy',ph_range);
        end
        plot_amp_phase(diff_c,[slices_diff_dirs],[-pi/2 pi/2],asp,params_in);diff_c=[];end
    close all
    
    %histogram of amp values
    params_sh.prefix='-A-';
    %params.colors=['r','b'];
    params_sh.xaxis='Amplitude Value';
    if shell_hist == 1
        if params_in.use_shells_from_first == 1,shell_sup=vol_series(:,:,:,1);disp('Using intial shells for remainder....');
        else shell_sup=vol_series(:,:,:,qq);disp('Calculating new shells....');
        end
        params_sh=density_shells_histogram(abs(crystal),shell_sup,nshells,[save_dir,'/',string_prefix,char(numbers_string(qq)),'/'],params_sh);
        hparams.amp_shell_vals{qq}=params_sh.shell_vals;
        params_sh.shell_vals=[];
    end

    %hist of phase values
    params_sh.prefix='-P-';
    %params.colors=['c','b'];        %used for the histogram colors
    params_sh.xaxis='Phase (2 \pi)';
    if shell_hist == 1
        if params_in.use_shells_from_first == 1,shell_sup=vol_series(:,:,:,1);disp('Using intial shells for remainder....');
        else shell_sup=vol_series(:,:,:,qq);disp('Calculating new shells....');
        end
        params_sh=density_shells_histogram((ph+pi)/2/pi,shell_sup,nshells,[save_dir,'/',string_prefix,char(numbers_string(qq)),'/'],params_sh);
        hparams.ph_shell_vals{qq}=params_sh.shell_vals;
        params_sh.shell_vals=[];
    end
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

%get cc fomr shells
if shell_hist == 1
    try
        [ hparams ] = calc_cc_shells(hparams);
    catch
        disp('<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>')
        disp('Error calculating cc for histogram shells....')
        disp('Try setting params.use_shells_from_first=1 ....')
    end
end
%save aligned arrays with the support
if params_in.save_aligned == 1
   disp('Saving aligned reconstructions....')
   aligned_arr=[save_dir,'/Aligned/'];
   make_directory({aligned_arr});
   save([aligned_arr,'temp_series.mat'],'temp_series')
   save([aligned_arr,'vol_series.mat'],'vol_series')
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
params.strain_sig=strain_sig;
params.skew_ph=skew_ph;
params.array_volume=tot_n;
params.array_size=size(crystal);
params.dir=top_dir;
params.surface_area=surface_area;

params.amp_std=std_dev;
params.amp_mean=amp_mean;

save([save_dir,'/','PARAMS.mat'],'params');

% plot some of the stuff and save images
lw=1.5;
font_size=25;
if plot_cold == 1,stuff.n_temps_cold=n_temps_cold;else stuff.n_temps_cold=0;end

stuff.lw=lw;
stuff.font_size=font_size;
stuff.xlabel=cxlabel;

if params_in.output_data_to_csv == 1,csvwrite([save_dir,'/xaxis.csv'],numbers);end

stuff.ylabel='Volume (\Delta V)';
stuff.save_name=[save_dir,'/Volume-T'];
stuff.color='blue';
plot_generic(numbers,volume,stuff)
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],volume);end

stuff.ylabel='Surface area (\Delta S)';
stuff.save_name=[save_dir,'/SurfA-T'];
stuff.color='blue';
plot_generic(numbers,surface_area/surface_area(1),stuff)
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],surface_area);end

if shell_hist == 1
    
    stuff.ylabel='CC';
    stuff.color='blue';
    
    try
        for qq=1:hparams.nshells
            stuff.save_name=[save_dir,'/CC-Amp-Shell-',num2str(qq)];    
            plot_generic(numbers,squeeze(hparams.cc_amp(:,qq)),stuff)  
            if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],hparams.cc_amp(:,qq));end
            stuff.save_name=[save_dir,'/CC-Ph-Shell-',num2str(qq)];    
            plot_generic(numbers,squeeze(hparams.cc_ph(:,qq)),stuff)
            if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],hparams.cc_ph(:,qq));end
            
        end    
    catch
        
        disp('<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>')
        disp('Error displaying cc for histogram shells....')
        disp('Try setting params.use_shells_from_first=1 ....')
    end
end

stuff.ylabel='Mean amplitude (\Delta \langle A \rangle)';
stuff.save_name=[save_dir,'/Amp-mean-T'];
stuff.color='blue';
plot_generic(numbers,amp_mean/amp_mean(1),stuff)
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],amp_mean);end
            

stuff.ylabel='Std amplitude (\Delta \sigma_{A}';
stuff.save_name=[save_dir,'/Amp-std-T'];
stuff.color='blue';
plot_generic(numbers,std_dev/std_dev(1),stuff)
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],std_dev);end
            

stuff.ylabel='Sphericity';
stuff.save_name=[save_dir,'/Sph-T'];
stuff.color='blue';
plot_generic(numbers,sphericity,stuff)
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],sphericity);end
            

stuff.ylabel='CC_{r}';
stuff.save_name=[save_dir,'/nxcc-T'];
stuff.color='red';
plot_generic(numbers,nxcc,stuff)
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],nxcc);end

stuff.ylabel='CC_{d}';
stuff.save_name=[save_dir,'/dxcc-T'];
stuff.color='green';
if no_data == 0, plot_generic(numbers,dxcc,stuff);
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],dxcc);end;end

%stuff.ylabel='Amplitude Spread (\sigma / \mu)';
%stuff.save_name=[save_dir,'/amp-spread'];
%stuff.color='b';
%plot_generic(numbers,std_dev./amp_mean,stuff)

stuff.ylabel='Phase skewness';
stuff.save_name=[save_dir,'/ph-skew'];
stuff.color='b';
plot_generic(numbers,skew_ph,stuff)
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],skew_ph);end

%stuff.ylabel='Phase (Avg. per pixel)';
%stuff.save_name=[save_dir,'/strain-T'];
%stuff.color='r';
%plot_generic(numbers,strain,stuff)

stuff.ylabel='Phase (radians)';
stuff.save_name=[save_dir,'/strain-1sig-T'];
stuff.color='b';
plot_generic(numbers,strain_sig,stuff)
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],strain_sig);end

stuff.ylabel='Phase (radians)';
stuff.save_name=[save_dir,'/strain-rms-T'];
stuff.color='b';
plot_generic(numbers,strain,stuff)
if params_in.output_data_to_csv == 1,csvwrite([stuff.save_name,'.csv'],strain);end

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

name=[save_dir,'/LO-Ph-X'];x_ph_lines_o=x_ph_lines;
if offset == 1,
    name=[name,'-OS'];
    x_ph_lines=add_offset_line_series(x_ph_lines,numel(amp_fs),.5,1);
    ph_max=max(x_ph_lines(:));ph_min=min(x_ph_lines(:));
end
axis_r=[xls/2-xls/4 xls/2+xls/4 1.1*ph_min 1.1*ph_max]; 
legend_str=numbers_string;
plot_amp_mult_color(x_ph_lines,1*size(x_ph_lines,2),lw,font_size,name,axis_r,legend_str,'Phase (radians)',[],[],stuff)

name=[save_dir,'/LO-Ph-Y'];y_ph_lines_o=y_ph_lines;
if offset == 1,
    name=[name,'-OS'];
    y_ph_lines=add_offset_line_series(y_ph_lines,numel(amp_fs),.5,1);
    ph_max=max(y_ph_lines(:));ph_min=min(y_ph_lines(:));
end
axis_r=[yls/2-yls/4 yls/2+yls/4 1.1*ph_min 1.1*ph_max]; 
legend_str=numbers_string;
plot_amp_mult_color(y_ph_lines,1*size(y_ph_lines,2),lw,font_size,name,axis_r,[],'Phase (radians)',[],[],stuff)

name=[save_dir,'/LO-Ph-Z'];z_ph_lines_o=z_ph_lines;
if offset == 1,
    name=[name,'-OS'];
    z_ph_lines=add_offset_line_series(z_ph_lines,numel(amp_fs),.6,1);
    ph_max=max(z_ph_lines(:));ph_min=min(z_ph_lines(:));
end
axis_r=[zls/2-zls/4 zls/2+zls/4 1.1*ph_min 1.1*ph_max]; 
legend_str=numbers_string;
plot_amp_mult_color(z_ph_lines,1*size(z_ph_lines,2),lw,font_size,name,axis_r,[],'Phase (radians)',[],[],stuff)

%plot a subset of the plots according to 
if isempty(plot_subset) ~= 1
    if isdir([save_dir,'/SubLO/']) ==0,mkdir([save_dir,'/SubLO/']);end
    
    %<< X, no offset
    name=[save_dir,'/SubLO/Sub-LO-Ph-X'];
    temp=x_ph_lines(:,plot_subset);
    ph_max=max(temp(:));ph_min=min(temp(:));
    axis_r=[xls/2-xls/4 xls/2+xls/4 1.1*ph_min 1.1*ph_max]; 
    legend_str=numbers_string(plot_subset);
    plot_amp_mult_color(x_ph_lines(:,plot_subset),1*size(x_ph_lines(:,plot_subset),2),lw,font_size,[name,'-OS'],axis_r,legend_str,'Phase (radians)',[],[],stuff)

    temp=x_ph_lines_o(:,plot_subset);
    ph_max=max(temp(:));ph_min=min(temp(:));
    axis_r=[xls/2-xls/4 xls/2+xls/4 1.1*ph_min 1.1*ph_max]; 
    plot_amp_mult_color(temp,1*size(temp,2),lw,font_size,name,axis_r,legend_str,'Phase (radians)',[],[],stuff);
    
    %<< Y  
    name=[save_dir,'/SubLO/Sub-LO-Ph-Y'];
    temp=y_ph_lines(:,plot_subset);
    ph_max=max(temp(:));ph_min=min(temp(:));
    axis_r=[yls/2-yls/4 yls/2+yls/4 1.1*ph_min 1.1*ph_max];    %1.1*ph_min 1.1*ph_max]; 
    plot_amp_mult_color(y_ph_lines(:,plot_subset),1*size(y_ph_lines(:,plot_subset),2),lw,font_size,[name,'-OS'],axis_r,legend_str,'Phase (radians)',[],[],stuff)

    temp=y_ph_lines_o(:,plot_subset);
    ph_max=max(temp(:));ph_min=min(temp(:));
    axis_r=[xls/2-xls/4 xls/2+xls/4  -1 2];
    plot_amp_mult_color(temp,1*size(temp,2),lw,font_size,name,axis_r,legend_str,'Phase (radians)',[],[],stuff);
    
    %<< Z
    name=[save_dir,'/SubLO/Sub-LO-Ph-Z'];
    temp=z_ph_lines(:,plot_subset);
    ph_max=max(temp(:));ph_min=min(temp(:));
    axis_r=[xls/2-xls/4 xls/2+xls/4 1.1*ph_min 1.1*ph_max]; 
    plot_amp_mult_color(z_ph_lines(:,plot_subset),1*size(z_ph_lines(:,plot_subset),2),lw,font_size,[name,'-OS'],axis_r,legend_str,'Phase (radians)',[],[],stuff)

    temp=z_ph_lines_o(:,plot_subset);
    ph_max=max(temp(:));ph_min=min(temp(:));
    axis_r=[zls/2-zls/4 zls/2+zls/4 1.1*ph_min 1.1*ph_max]; 
    plot_amp_mult_color(temp,1*size(temp,2),lw,font_size,name,axis_r,legend_str,'Phase (radians)',[],[],stuff);
    
    
end


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
%jclark

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
    %saveas(fh, save_name, 'epsc');
    %print(fh, '-dpng','-r600', save_name);
    [save_dir save_name]=extract_dir_from_string(save_name);
    save_figure_mult(fh,save_dir,save_name,1)
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


%saveas(fh, name, 'epsc');
%print(fh, '-dpng','-r600', name);
[save_dir save_name]=extract_dir_from_string(name);
save_figure_mult(fh,save_dir,save_name,1)

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

function sz = plot_amp_phase(pn,save_dir,ph_range,sz,params_in)

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

if params_in.interp_slices == 1
    imagesc(imresize(ph(yy(1):yy(2),xx(1):xx(2)),2));
else
    imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
end

h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
axis equal
%caxis(ph_range)
%saveas(fh, [fdir,'Amp-xy'],'epsc');
%print(fh, '-dpng','-r300', [fdir,'Amp-xy']);
[xdir xname]=extract_dir_from_string(fdir);
save_figure_mult(fh,xdir,[xname,'Amp-xy'],1)

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(amp,'xz' );

if params_in.interp_slices == 1
    imagesc(imresize(ph(zz(1):zz(2),xx(1):xx(2)),2));
else
    imagesc(ph(zz(1):zz(2),xx(1):xx(2)));
end

h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
axis equal
%caxis(ph_range)
%saveas(fh, [sdir,'Amp-xz'],'epsc');
%print(fh, '-dpng','-r300', [sdir,'Amp-xz']);
[xdir xname]=extract_dir_from_string(sdir);
save_figure_mult(fh,xdir,[xname,'Amp-xz'],1)


fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(amp,'zy' );

if params_in.interp_slices == 1
    imagesc(imresize(ph(yy(1):yy(2),zz(1):zz(2)),2));
else
    imagesc(ph(yy(1):yy(2),zz(1):zz(2)));    
end

h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
axis equal
%caxis(ph_range)
%saveas(fh, [tdir,'Amp-zy'],'epsc');
%print(fh, '-dpng','-r300', [tdir,'Amp-zy']);
[xdir xname]=extract_dir_from_string(tdir);
save_figure_mult(fh,xdir,[xname,'Amp-zy'],1)


if isreal(pn) ~= 1
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    %ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
    ph = extract_3D_slice(phase,'xy' );
    
    if params_in.interp_slices == 1
        imagesc(imresize(ph(yy(1):yy(2),xx(1):xx(2)),2));
    else
        imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
    end
    
    h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
    axis equal
    caxis(ph_range);
    %saveas(fh, [fdir,'Ph-xy'],'epsc');
    %print(fh, '-dpng','-r300', [fdir,'Ph-xy']);
    [xdir xname]=extract_dir_from_string(fdir);
    save_figure_mult(fh,xdir,[xname,'Ph-xy'],1)

    
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    ph = extract_3D_slice(phase,'xz' );
    
    
    if params_in.interp_slices == 1
        imagesc(imresize(ph(zz(1):zz(2),xx(1):xx(2)),2));
    else
        imagesc(ph(zz(1):zz(2),xx(1):xx(2)));
    end
    
    h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
    axis equal
    caxis(ph_range);
    %saveas(fh, [sdir,'Ph-xz'],'epsc');
    %print(fh, '-dpng','-r300', [sdir,'Ph-xz']);
    [xdir xname]=extract_dir_from_string(sdir);
    save_figure_mult(fh,xdir,[xname,'Ph-xz'],1)

    
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    ph = extract_3D_slice(phase,'zy' );
    
    if params_in.interp_slices == 1
        imagesc(imresize(ph(yy(1):yy(2),zz(1):zz(2)),2));
    else
        imagesc(ph(yy(1):yy(2),zz(1):zz(2)));
    end
    
    h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
    axis equal
    caxis(ph_range);
    %saveas(fh, [tdir,'Ph-zy'],'epsc');
    %print(fh, '-dpng','-r300', [tdir,'Ph-zy']);
    [xdir xname]=extract_dir_from_string(tdir);
    save_figure_mult(fh,xdir,[xname,'Ph-zy'],1)

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

function params=set_params_defaults(params)
try
    params.no_data;
catch
    params.no_data=1;    %don't do the data, its takes time, use data = 0
end

try
    params.avgr;
catch
    params.avgr=3;                           %region to average over, nxn
end

try
    params.offset;
catch
    params.offset=1;                         %plot the lineouts with an offset?, yes=1
end
    
try
    params.line_out_cent_offset;
catch
    params.line_out_cent_offset=[0,0,0];     %offset of where to take the lineouts from
end

try
    params.shell_hist;
catch
    params.shell_hist=0;             %do the shell histograms, yes=1
end

try
    params.amp_norm;
catch
    params.amp_norm=0;               %normalise amplitudes of reconstruction to that of the first
end                            %such that max(first)=1, and then %sum(second)=sum(first).....
                                 
%<<< output params >>>%
try
    params.conj_first;
catch
    params.conj_first=0;            %conjugate the reference? (=1 yes, 0=no)
end

try
    params.lab;
catch
    params.lab=1;                   %load the transformed reconstruction
end

try
    params.vtk_out;
catch
    params.vtk_out=1;               %output vtk file (will be correctly orientated)
end

try
    params.rem_ramp;
catch
    params.rem_ramp=1;              %remove ramp in sample plane, =1 yes
end

try
    params.upsize;
catch
    params.upsize=3;                %upsample factor
end

try
    params.z_phase;
catch
    params.z_phase=1;               %set mean to zero (1=yes)
end

try
    params.apply_sup;
catch
    params.apply_sup=0;             %multiply by support (vtk output)
end

try
    params.apply_sup_slices;
catch
    params.apply_sup_slices=0;      %apply sup for slices as well
end

try
    params.output_ind_slices;
catch
    params.output_ind_slices=0;
end

try
    params.output_ind_N;
catch
    params.output_ind_N=32;
end

%<<< cold params >>>
try
    params.plot_cold;
catch
    params.plot_cold=0;                %set =1 to plot the cooled down temperature 30 if present
end

try
    params.cold_string;
catch
    params.cold_string='a/r';          %s string that is in the file name of the cold temps that distniguishes them from the rest
end

%<<< Misc params >>>>
try
    params.sort_numerically;
catch
    params.sort_numerically=0;       %will sort according to number.  do not use if heating series has coold down temps, will lose the order
end

try
    params.nnumbs_out;
catch
    params.nnumbs_out=[];            %number of digits in the numbers for save name output, use [] to turn off
end

%<<< Plots params >>>>
try
    params.auto_size_plot;
catch
    params.auto_size_plot=0;        %auto size (=0), use all =1
end
    
try
    params.hist_range;
catch
    params.hist_range=.5;           %fraction of 2 pi
end

try
    params.ph_bins;
catch
    params.ph_bins=100;
end

try
    params.ph_range;
catch
    params.ph_range=[-pi/2,pi/2];           %range for output image
end

try
    params.offset_amp;
catch
    params.offset_amp=.65;              %offset to plot amp histogam as fraction of max
end

try
    params.offset_ph;
catch
    params.offset_ph=.6;                %offset to plot ph histogam as fraction of max
end

try
    params.vsig;
catch
    params.vsig=0.5;                    %calc values for the volume
end
    
try
    params.vth;
catch
    params.vth=0.1;
end

try
    params.stuff.bold_plot;
catch
    params.stuff.bold_plot=0;
end

try
    params.stuff.dots_only;
catch
    params.stuff.dots_only=0;
end

try
    params.custom_xaxis;
catch
    params.custom_xaxis=[];          %use a custom axis, leave as [] for nothing
end
   
try
    params.cxlabel;
catch
    params.cxlabel='Temperature (C)';%;
end

try
    params.points;
catch
    params.points=[30,47,36];  %x,y,z, J9
end

try
    params.calc_diff;
catch
    params.calc_diff=0;
end

try
    params.calc_diff_slices;
catch
    params.calc_diff_slices=0;
end

try
    params.interp_slices;
catch
    params.interp_slices=0;
end

try
    params.upsamp_vtk;
catch
    params.upsamp_vtk=1;
end


try
    params.box_vtk;
catch
    params.box_vtk=0;
end

try
    params.save_aligned;
catch
    params.save_aligned=0;
end

try
    params.save_aligned_name;
catch
    params.save_aligned_name='';
end

try
    params.nshells;
catch
    params.nshells=2;
end

try
    params.shell_maxy;
catch
    params.shell_maxy=1000;
end

try
    params.use_shells_from_first;
catch
    params.use_shells_from_first=0;
end

try
    params.output_shell_to_csv;
catch
    params.output_shell_to_csv=0;
end

try
    params.output_data_to_csv;
catch
    params.output_data_to_csv=0;
end

try
    params.crop_vtk;
catch
    params.crop_vtk=[];
end

end
