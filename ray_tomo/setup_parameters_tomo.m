% Setup_parameters for raytomo
%
% https://github.com/jbrussell
addpath('../');
addpath('../functions');
addpath('./tomo_functions');
% parameters.qpath = './MINEOS_qfile/';

parameters.station_list = '../../sta_locs.txt';

% Lat Lon
parameters.lalim = [32.0 37.0] ;
parameters.lolim = [-77.0 -71.0];
parameters.gridsize = 0.25; % degrees?
parameters.agebins = [165:5:175];
parameters.bathybins = [-9000 :5000: 1000];
parameters.gridsize_azi = 0.25; %3; %1.5; % gridsize for 2D azimuthal anisotropy (degrees)
parameters.r = 0.03; %0.01; % controls color bar [avgv(1-r) avgv(1+r)]

% Smoothing parameters
parameters.smweight0 = 100; % isotropic second derivative smoothing
parameters.smweight0_azi = 1e3; %1000; % anisotropic second derivative smoothing
parameters.flweight0_azi = 1000; %1000; % anisotropic first derivative flatness

% parameters for the tomography (QC)
parameters.raydensetol=deg2km(parameters.gridsize)*0.25; %deg2km(parameters.gridsize); %deg2km(parameters.gridsize)*2;
parameters.raydensetol_azi=deg2km(parameters.gridsize_azi)*0.25; %deg2km(parameters.gridsize)*2;
parameters.fiterrtol = 2;   % error allowed in the wavelet fitting
parameters.dterrtol = 4;    % largest variance of the inversion error allowed
parameters.maxerrweight = 5; % Maximum error weight
parameters.polyfit_dt_err = 2; % (s) dt error greater than this, weighted 0