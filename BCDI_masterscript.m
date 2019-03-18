% in this script the scans are read and the parameters are initiated

global  X Y Z d2_bragg ki_o kf_o

addpath(genpath('/Users/ialmazn/Box Sync/forIrene/openspec-1.4'));
addpath(genpath('/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Analysis_end_of_beamtime'));
addpath(genpath('/Users/ialmazn/Documents/MATLAB/ptycho/m_scripts/'));
addpath(genpath('/Users/ialmazn/Box Sync/forIrene/Matlab-Phasing-27May-2014'));
% read spec file:

specnum = 92;

[specscan, errors] = openspec('Stephenson316a.spec',specnum );

% the geometry of the experiment:
BCDI_diff_vectors;


% load the tif files into a matrix:

filename =  ['Stephenson316a_S' num2str(specnum,'%04d') '_'];
directname =  ['Stephenson316a_S' num2str(specnum,'%04d') '/'];

mindata = 3; % sets a threshold for background 
%[data_exp,rock,imgs] = BCDI_read_center_scans(filename,directname,specscan);
[data_exp,rock,imgs] = BCDI_read_center_pad_scans(filename,directname,specscan,mindata);

arrysize = size(imgs,1);
arrysize3 = size(imgs,3);


%probe = ones(size(X));

%BCDI_ph_retrieval;

%onthefly_GA;

onthefly_GA_Irene;


