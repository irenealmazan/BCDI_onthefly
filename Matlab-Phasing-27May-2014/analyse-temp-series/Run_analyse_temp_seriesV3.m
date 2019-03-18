function Run_analyse_temp_seriesV2
%jclark
%wrapper for running the analyse temp series function
%duplicate this file and use it to set the params for different series

top_dir='/Users/jesseclark/Documents/MATLAB/data_analysis/Au0812/J2/C6/';   %the top dir to search below, leave =[] to use the path of this file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.prefix='Rec-A-151*-GAHIO-AVG-5-sS-8-5-200-CVl-SW';
params.string_prefix='T';
params.outdir='TempSeries';

%set custom_file_names=[] to use the automtaic version. don't include the
%top dir in the names
params.custom_file_names=make_string_from_recursions('Rec-',{'1092','1107'},'-100-GAHIO-AVG-4-sS-16-10-200-NM-SW');%'1187','1170','1222','1291','1354','1337','1375','1431','1447','1478','1490'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.no_data=1;    %don't do the data, its takes time, use data = 0

%<<< line out params >>>%
params.avgr=3;                           %region to average over, nxn
params.offset=1;                         %plot the lineouts with an offset?, yes=1
params.line_out_cent_offset=[0,0,0];     %offset of where to take the lineouts from

params.shell_hist=0;             %do the shell histograms, yes=1
params.nshells=4;                %number of shells to use
params.shell_maxy=5000;

params.use_shells_from_first=1;  %use shells defind from first member
params.output_shell_to_csv=1;    %save each shell to csv 

params.output_data_to_csv=0;     %output all data to csv

params.amp_norm=0;               %normalise amplitudes of reconstruction to that of the first
                                 %such that max(first)=1, and then %sum(second)=sum(first).....
                                 
%<<< output params >>>%
params.conj_first=1;            %conjugate the reference? (=1 yes, 0=no)
params.lab=1;                   %load the transformed reconstruction
params.vtk_out=1;               %output vtk file (will be correctly orientated)
params.rem_ramp=1;              %remove ramp in sample plane, =1 yes
params.upsize=3;                %upsample factor
params.z_phase=1;               %set mean to zero (1=yes)
params.apply_sup=1;             %multiply by support (vtk output)

params.output_ind_slices=1;    %output slices?
params.output_ind_N=16;             %number of slices to output

params.apply_sup_slices=1;      %apply sup for slices as well

params.vtk_markers=1;           %put markers in the corners to make all arrays display the same

params.vtk_zero2nan=1;          %save zeros as Nans
%<<< cold params >>>
params.plot_cold=1;                %set =1 to plot the cooled down temperature 30 if present
params.cold_string='a/r';          %s string that is in the file name of the cold temps that distniguishes them from the rest

%<<< Misc params >>>>
params.sort_numerically=0;       %will sort according to number.  do not use if heating series has coold down temps, will lose the order
params.nnumbs_out=[];            %number of digits in the numbers for save name output, use [] to turn off

%<<< Plots params >>>>
params.auto_size_plot=0;        %auto size (=0), use all =1
params.hist_range=.5;           %fraction of 2 pi
params.ph_bins=100;

params.ph_range=[-pi/2,pi/2];           %range for output image

params.offset_amp=.65;              %offset to plot amp histogam as fraction of max
params.offset_ph=.6;                %offset to plot ph histogam as fraction of max

params.vsig=1.0;                    %calc values for the volume
params.vth=0.1;

params.stuff.bold_plot=0;
params.stuff.dots_only=0;

params.custom_xaxis=[];          %use a custom axis, leave as [] for nothing, should be a number
params.cxlabel='Temperature (C)';%;

params.points=[30,47,36];  %x,y,z, J9

params.calc_diff=0;

params.interp_slices=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(top_dir)
    name_of_this_file='Run_analyse_temp_seriesV3';
    dir_file=which(name_of_this_file);    %finds the current folder for the phasing
    top_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
end
    
analyse_temp_seriesV3(top_dir,params.prefix,params);

make_movies_T_analysis([top_dir,params.outdir,'/'],5,[])

end
