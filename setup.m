addpath(genpath([pwd,'\bestfitparams']))
addpath(genpath([pwd,'\data']))
addpath(genpath([pwd,'\model']))
addpath(genpath([pwd,'\visuals']))
addpath(genpath([pwd,'\zIP']))

% Set figure properties
set(0, 'DefaultLineLineWidth', 2);
set(groot,'DefaultAxesFontSize',16);
set(0,'Defaultfigurecolor',[1 1 1]);

% Constant parameters used throughout
data_filenames = {'dosSantos2017_old','dosSantos2017_young'};
data_types = {'Older','Younger'};
colors = [0    0.447    0.741
           0.850    0.325    0.098];
alpha = 1e6;

% Simulation parameters
input_struct.simFreq_Hz         = 1000;
input_struct.sampFreq_Hz        = 1000;
input_struct.simDuration_s      = 1000;
input_struct.motorNoiseLvL_Nm   = 1;
input_struct.noise_type         = 'w';
