% Setup_parameters for ambient noise processing
%
% NJA, 4/2/2016
% JBR, 6/16/2016

addpath('./functions/');
addpath('./functions/calc_Rayleigh_disp/');

%%% --- Paths to important files --- %%%
parameters.workingdir = pwd;
parameters.workingdir = [parameters.workingdir,'/'];

parameters.datapath = ['/data/irma6/jrussel/YoungPacificORCA/SAC_1Hz/']; %'../nomelt_data_5sta/';
% parameters.datapath = ['/data/irma6/jrussel/YoungPacificORCA/SAC_1Hz_Zcorr_tiltcomp/']; %'../nomelt_data_5sta/';

parameters.PZpath = '../INSTRUMENT/';
parameters.ccfpath = './ccf/'; %['/Volumes/LaCie/AmbNoise_allCCFs_nomelt_molly/'];
parameters.figpath = [parameters.workingdir,'figs/'];
parameters.seis_path = [parameters.workingdir,'seismograms/'];
%parameters.orientation_path = '../OBS_Orientation/Orientation_Noise/';
parameters.orientation_path = '/Users/jrussel/RESEARCH/PROJ_YoungPacificORCA/DATA/ORCA_orientations.txt';
% parameters.stalist = textread([parameters.workingdir,'stalist'],'%s\n');
[stalist, stalat, stalon, staz] = textread(['stalist_good12.txt'],'%s %f %f %f\n');
parameters.stalist = stalist;
parameters.stalat = stalat;
parameters.stalon = stalon;
parameters.staz = staz;
parameters.nsta = length(parameters.stalist);

%%% --- Parameters to build up gaussian filters --- %%%
parameters.min_width = 0.18; %0.06
parameters.max_width = 0.30; %0.1

%%% --- Parameters for initial processing --- %%%
parameters.dt = 1; % sample rate
parameters.comp = 'BH'; % component
parameters.mindist = 20; % min. distance in kilometers

%%% --- Parameters for ccf_ambnoise --- %%%
parameters.winlength = 3; %hours

%%% --- Parameters for fitbessel --- %%%
parameters.npts = parameters.winlength*3600;

%%% --- Parameters for OBS orientations --- %%%
% all_stalist = { 'B01', 'B02', 'B04', 'B05','B06','B08', 'B11', 'B13', 'B16','B17','B19','B22','B23', 'B24','B25','B26'  };
