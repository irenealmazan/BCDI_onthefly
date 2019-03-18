global d2_bragg X Y Z ki_o kf_o
warning off;
addpath(genpath('/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Analysis_end_of_beamtime'));
%addpath(genpath('..\..\..\m_scripts\'));

%Meganrun = 1;
%if Meganrun == 1
addpath(genpath('/Users/ialmazn/Documents/MATLAB/ptycho/m_scripts/'))
%addpath(genpath('C:\Users\moh9\Documents\NSLS2_March17\Data'))
%end


addCOMSOL = 0;
%experimental setup details
pixsize = 55; %microns, for Merlin
lam = etolambda(10400)*1e-4;
a_latparam = 4.24e-4; %for InGaAs, 14.3% Ga fraction
c_latparam = 6.93e-4; 
mplane_spacing = sind(60)*a_latparam;
aplane_spacing = a_latparam/2;
q_mplane = 2*pi/(mplane_spacing*1e4);
q_aplane = 2*pi/aplane_spacing *1e-4;
q_cplane = 2*pi/c_latparam * 1e-4;
%th_mplane110 = asind(lam/(2*mplane_spacing));
%th_mplane220 = asind(2*lam/(2*mplane_spacing));
%th_aplane112 = asind(lam/(2*a_latparam/2));
meshdata = 0;
%ZP details
zpdiam = 240;
outerzone = .04;
bsdiam = 50;
binaryprobe_flag =0; %use binary probe defined at FWHM of probe

%NW details
corewidth= .2; %edge-to-edge distance %NW 22
%corewidth = 2*cosd(50)* .087; %take estimate from SEM image
length = 2;
phi=0; %misalignment of NW from vertical

%scattering condition:
scatgeo = 2110; %for strain image
%scatgeo = 1010; %for stacking-fault image

%cutrad = .12;
cutrad = .06;
edgepad = 1.25;
usesimI = 1; 

if scatgeo==1010 %SF
    
    Npix = 516;
    detdist = 0.33e6;
    d2_bragg = detdist * lam /(Npix*pixsize);
    depth = 100; 
    defocus = 0;

    %expected values for m-plane
    th = -9.3;
    del = -18.6; %in plane
    gam = 0; %out of plane

    addNWstrain = 0;
    addNWsf = 0;
    
    
    %NW 20 (on Cr grid line), R1 th/area scan, p100
    fly2Dscanlist = [25304:4:25360]; fly2Danglist = [-9.42:.02:-9.13];    
    fly2Dscanlist = fly2Dscanlist(2);
    fly2Danglist = fly2Danglist(2);
    thscanlist = fly2Dscanlist; 
    thscanvals = fly2Danglist;
    thBragg = -9.4;
    cen_pix = [244 52]; %row colum
    %NW 20 (on Cr grid line), R2 th/area scan, p102
    %fly2Dscanlist = [0:14]*4+25389; fly2Danglist = [-9:.02:-8.72];
    %thscanlist = fly2Dscanlist; thscanvals = fly2Danglist;
    %thBragg = -8.9;
    
    %NW 22 , whole wire th/area scan, p106
    %fly2Dscanlist = [25573:2:25621]; fly2Danglist = [-9.96:.04:-9.0];
    %thscanlist = fly2Dscanlist; thscanvals = fly2Danglist;
    %thBragg = -9.7;

    thscanvals = thscanvals-thBragg;
    
    noiseflag = 1;
    mncntrate = 1;
    useflyscan = 1;
        
    edgepad = 1.3; %for support during phase ret.
    %display('Predefined Edgepad temp')
end

if scatgeo==2110 %strain 

    %most common scan conditions - used unless otherwised changed
    lineup = 0;  %no 2D line up scan performed
    %det vals for, a-plane peak
    th = 73.3;
    del = -32.6; %in plane
    gam = 0; %out of plane


    
    
    %NW 20 - fine 2D at many angles
    fly2Dscanlist=[(24938+((1:8)-1)*3) (24975+((1:5)-1)*3) (24995+((1:25)-1)*3)];
    fly2Danglist =  [72.8:1.0/64:73.8-1.0/64];%[73.0:0.01:73.74];%[(73.0+((1:8)-1)*0.02) (73.16+((1:5)-1)*0.02) (73.26+((1:25)-1)*0.02)];
    %thscanlist = fly2Dscanlist(1:1:numel(fly2Danglist));
    %thscanvals = fly2Danglist(1:1:numel(fly2Danglist));
    
    % side region
    %{
    thscanlist = fly2Dscanlist(24:1:28);
    thscanvals = fly2Danglist(24:1:28);
    thBragg = 73.5;
    cen_pix = [135,155];
    %}
    %Center region
    %%{
    thscanlist = fly2Danglist;%fly2Dscanlist(15:17);%fly2Dscanlist(15);%
    thscanvals = fly2Danglist;%fly2Danglist(15:17);%fly2Danglist(15);%
    thBragg = 73.3; 
    cen_pix = [89 150]; %x,y, unrotated as-read image
    %cen_pix = [(217+1)/2, (340+1)/2];
    %}
    
 
   display('thscanvals, thBragg, cen_pix, defined by outside script')
    thscanvals = thscanvals-thBragg;

    Npix = 256;
    detdist = 0.529e6; %um - using scan#: 24915, 24918 w/ vert change of 0.5 deg. 
    d2_bragg = detdist * lam /(Npix*pixsize);
    depth = size(thscanvals,2);%75;%60;%60;
    defocus = 0;
    
    %{
    % calculate the qbragg vector:
    kf = [0 0 1]; % detector frame
    ki = [0 0 1]; % lab frame
    kmag = 2*pi/lam;
    
    Ry = [cosd(-del) 0 sind(-del);
        0 1 0;
        -sind(-del) 0 cosd(-del)];
    
    ki = (Ry * ki.').';
    
    Rx=[1 0 0;
        0 cosd(gam) -sind(gam);
        0 sind(gam) cosd(gam)];
    ki = (Rx * ki.').'; % detector frame
    
    qbragg = kf-ki; % detector frame
    
    ki_o = ki;
    kf_o = kf;
    
    % calculate the resolution in the theta direction
    dth_ini = fly2Danglist(1) - thBragg;
    
    
    Ry = [cosd(-dth_ini) 0 sind(-dth_ini);
        0 1 0;
        -sind(-dth_ini) 0 cosd(-dth_ini)];
    
    ki = (Ry * ki_o.').';
    kf = (Ry * kf_o.').';
    
    dq_ini(1,:) = (kf-ki)-qbragg;
    
    dth_fin = fly2Danglist(end) - thBragg;
    
    
    Ry = [cosd(-dth_fin) 0 sind(-dth_fin);
        0 1 0;
        -sind(-dth_fin) 0 cosd(-dth_fin)];
    
    ki = (Ry * ki_o.').';
    kf = (Ry * kf_o.').';
    
    dq_fin = (kf-ki)-qbragg;
    
    qz_pixel_size = abs((dq_fin(3) - dq_ini(3))/numel(fly2Danglist));

    th_pixel_size = 2*pi/abs((dq_fin(3) - dq_ini(3)))*1e-3;
    %}
    
    addNWstrain = 0;
    addNWsf = 0;

    mncntrate = 1;
    %usesimI = 1;
    if(usesimI) display('USING SIMULATED DATA');end
end

%%
NW_diff_vectors_BCDI_v2; % does both the vectors ki and kf and creates the object
%NW_diff_vectors_BCDI;
%NW_make_InGaAs_nocoreshell_BCDI;
probe = ones(size(X));
%NW_diff_plot_vectors_BCDI;
%NW_diff_vectors_Irene;




%NW_calc_rocking_curve;
%NW_calc_dp;
%NW_calc_dp_Irene;
number_angles_distort = 0; % number of angles to distort
NW_calc_dp_BCDI;

if(usesimI) 
    randomIniGuess = 0;
    NW_use_simdata; 
end
%NW_add_dp_noise;

%NW_theta_annealing_Irene(probe, NW,data_exp,dth_disp, thBragg);

%NW_theta_annealing_Irene_test_v2(probe, NW,data_exp,dth_disp, thBragg);

onthefly_GA_v2;

%test_FT_InvFT_Irene;

return;

NW_ph_retrieval_Irene_BCDI;
return;
NW_test_gradient;
return;
%NW_ph_retrieval_Irene;



return

NW_test_gradient_2DBPP;
NW_BPP_ph_retrieval;