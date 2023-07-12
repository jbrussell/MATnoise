% SETUP PARAMETERS FOR BEAMFORMING
addpath('./functions_beamforming/'); addpath('../functions/');


%% Data paths

% Input data
parameters.comp = {'ZZ'}; %'ZZ'; %'PP'; 'RR'; 'TT';
parameters.windir = 'window3hr';
% parameters.windir = 'window3hr_Zcorr_tiltcomp'; % Tilt & compliance corrected

parameters.ccfpath = '../ccf/';
parameters.workingdir = './';
parameters.figpath = [parameters.workingdir,'figs/'];
[stalist, stalat, stalon, staz] = textread(['../stations_oldorca.txt'],'%s %f %f %f\n');
parameters.stalist = stalist;
parameters.stalat = stalat;
parameters.stalon = stalon;
parameters.staz = staz;
parameters.nsta = length(parameters.stalist);

%% PARAMETERS
% Time window to consider in beamforming
parameters.t_min = -200;
parameters.t_max = 200;

% Periods to average over
parameters.per_min = 18; % [sec] minimum period
parameters.per_max = 25; % [sec] maximum period
parameters.Npers = 30; % number of periods to consider

% Slowness values to search over
parameters.s_min = 1/5; % s/km
parameters.s_max = 1/2.5; % s/km
parameters.Nslow = 100/2; % s/km

% Back-azimuth values to search over
parameters.Nbaz = 360/2;

