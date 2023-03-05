% Setup_parameters for ambient noise processing
%
% NJA, 4/2/2016
% JBR, 6/16/2016

addpath('./functions/');
addpath('./functions/calc_Rayleigh_disp/');

%%% --- Paths to important files --- %%%
parameters.workingdir = [pwd,'/'];

parameters.datapath = ['/data/irma6/jrussel/YoungPacificORCA/SAC_1Hz/']; %'../nomelt_data_5sta/';

parameters.PZpath = '../INSTRUMENT/'; % path to RESP files containing poles and zeros
parameters.ccfpath = './ccf/';
parameters.figpath = [parameters.workingdir,'figs/'];
parameters.seis_path = [parameters.workingdir,'seismograms/'];
parameters.orientation_path = './OBS_orientations.txt'; % Column 1: station name;   Column 2: H1 degrees CW from N
[stalist, stalat, stalon, staz] = textread(['stalist_good12.txt'],'%s %f %f %f\n');
parameters.stalist = stalist;
parameters.stalat = stalat;
parameters.stalon = stalon;
parameters.staz = staz;
parameters.nsta = length(parameters.stalist);

%%% --- Parameters to build up gaussian filters --- %%%
parameters.min_width = 0.18;
parameters.max_width = 0.30;

%%% --- Parameters for initial processing --- %%%
parameters.dt = 1; % sample rate
parameters.comp = 'BH'; % component
parameters.mindist = 20; % min. distance in kilometers

%%% --- Parameters for ccf_ambnoise --- %%%
parameters.winlength = 3; %hours
parameters.Nstart_sec = 50; % number of sections to offset start of seismogram

%%% --- Parameters for fitbessel --- %%%
parameters.npts = parameters.winlength*3600;

%%% --- Parameters for using Radon Transform picks --- %%%
parameters.path_LRT_picks = './mat-LRTdisp/LRT_picks/';
