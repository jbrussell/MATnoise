% Setup_parameters for ambient noise processing
%
% NJA, 4/2/2016
% JBR, 6/16/2016

addpath('./functions/');
addpath('./functions/calc_Rayleigh_disp/');

%%% --- Paths to important files --- %%%
parameters.workingdir = pwd;
parameters.workingdir = [parameters.workingdir,'/'];

parameters.datapath = ['/Users/jrussel/RESEARCH/PROJ_ENAM/DATA/NOISE/YO_LH/']; %'../nomelt_data_5sta/';

parameters.PZpath = '../INSTRUMENT/';
parameters.ccfpath = './ccf/';
parameters.figpath = [parameters.workingdir,'figs/'];
parameters.seis_path = [parameters.workingdir,'seismograms/'];
parameters.orientation_path = '/Users/jrussel/RESEARCH/PROJ_ENAM/DATA/ORIENTATIONS/YO_orientations.txt';
[stalist, stalat, stalon, staz] = textread(['sta_locs.txt'],'%s %f %f %f\n');
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
parameters.comp = 'HH'; % component
parameters.mindist = 20; % min. distance in kilometers

%%% --- Parameters for ccf_ambnoise --- %%%
parameters.winlength = 3; %hours

%%% --- Parameters for fitbessel --- %%%
parameters.npts = parameters.winlength*3600;
