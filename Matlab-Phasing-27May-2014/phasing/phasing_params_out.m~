function phasing_params_out( params ,dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


deg2rad = pi/180.0000;
%# Detector and scan variables
binx=params.binning(1);
biny=params.binning(2);
det_px=params.det_px;

delta_deg=params.delta;
gam_deg=params.gam;

ang_n=params.dth;

pixelx=binx*(det_px);  %#Total pixel size including binning (in meters)
pixely=biny*(det_px);
arm = params.arm       ;    %#Detector distance from sample (units are m)
lam = params.lam    ;      %#Wavelength of x-ray (units in nm)
delta = delta_deg * deg2rad ;   %#Detector angles
gam =  gam_deg * deg2rad;
dth = ang_n * deg2rad;      %#step values for th translation
dtilt = params.dtilt *deg2rad;


% open a file for writing
fid = fopen([dir,'phasingparams.py'], 'w');

% print a title, followed by a blank line
fprintf(fid, 'import math \n');
fprintf(fid, '# Constants \n');
fprintf(fid, 'pi=math.pi \n');
fprintf(fid, 'deg2rad = pi/180.0000 \n');
fprintf(fid, '# Detector and scan variables \n');
fprintf(fid, 'det_px=');
fprintf(fid, '%20.15f\n', det_px);
fprintf(fid, 'binx=');
fprintf(fid, '%4.2f\n', binx);
fprintf(fid, 'biny=');
fprintf(fid, '%4.2f\n', biny);
fprintf(fid, 'delta_deg=');
fprintf(fid, '%6.3f\n', delta_deg);
fprintf(fid, 'gam_deg=');
fprintf(fid, '%6.3f\n', gam_deg);
fprintf(fid, 'pixelx=binx*det_px \n');
%fprintf(fid, '%20.15f\n', pixelx)
fprintf(fid, 'pixely=biny*det_px \n');
%fprintf(fid, '%20.15f\n', pixely)
fprintf(fid, 'arm=');
fprintf(fid, '%20.15f\n', arm);
fprintf(fid, 'lam=');
fprintf(fid, '%20.15f\n', lam);
fprintf(fid, 'delta=delta_deg*deg2rad \n');
fprintf(fid, 'gam=gam_deg*deg2rad \n');
fprintf(fid, 'dth=');
fprintf(fid, '%20.15f\n', dth);
fprintf(fid, 'dtilt=');
fprintf(fid, '%20.15f\n', dtilt);
% print values in column order
% two values appear on each row of the file
%fprintf(fid, '%f  %f\n', y);


%fprintf(fid, '\nparam_x=');
%pfprintf(fid, '%f  %f\n', x);


fclose(fid);


end

