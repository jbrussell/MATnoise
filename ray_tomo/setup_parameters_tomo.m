% Setup_parameters for raytomo
%
% https://github.com/jbrussell
addpath('../');
addpath('../functions');
addpath('./tomo_functions');
% parameters.qpath = './MINEOS_qfile/';

parameters.station_list = '../stalist_good12.txt';

% Lat Lon
parameters.lalim = [-9 -3] ;
parameters.lolim = [-136 -130];
parameters.gridsize = 0.25; % degrees?
parameters.agebins = [165:5:175];
parameters.bathybins = [-9000 :5000: 1000];
parameters.gridsize_azi = 0.5; %3; %1.5; % gridsize for 2D azimuthal anisotropy (degrees)
parameters.r = 0.03; %0.01; % controls color bar [avgv(1-r) avgv(1+r)]

% Parameters for reading XSP files
parameters.comp = {'ZZ'}; % component
parameters.xspdir = 'phv_dir'; % directory name
parameters.windir = 'window3hr'; % window type
parameters.N_wl = 1; % number of wavelengths to consider during bessel fit measurements
parameters.frange = [1/10 1/5]; % [Hz]
parameters.per_ind = [1:12]; % index of periods to consider
% QC parameters
parameters.snr_tol = 3; % minimum signal-to-noise
parameters.r_tol_min = 0; % [km] minimum station separation
parameters.r_tol_max = 600; % [km] maximum station separation
parameters.err_tol = inf; % maximum misfit of bessel fit between observed and synthetic
parameters.is_rtolmin_wavelength = 1; parameters.wl_fac = 1.0; % determine distance tolerance by wavelength?
parameters.dep_tol = [0 0]; % [sta1 sta2] OBS depth tolerance
parameters.is_raydensity_thresh = 0; % Apply raydensity threshold to wipe out poorly constrained grid cells?
parameters.min_dep = 9999; %-3500 for min station depth to use

% 1-D Anisotropy parameters (1Dazi scripts)
parameters.azi_bin_deg = 20; % (degrees) size of azimuthal data bin
% Norm damping for azimuthal anisotropy
parameters.damp_azi = [1 1 1e10 1e10]; % [2c 2s 4c 4s] % Damping individual parameters
parameters.aziweight = 1; % global weight

% Smoothing parameters (2D phv and anisotropy)
parameters.smweight0 = 100; % isotropic second derivative smoothing
parameters.smweight0_azi = 1e3; %1000; % anisotropic second derivative smoothing
parameters.flweight0_azi = 1000; %1000; % anisotropic first derivative flatness
parameters.damp0_azi = 0; % anisotropic norm damping

% parameters for the tomography (QC)
parameters.raydensetol=deg2km(parameters.gridsize)*0.25; %deg2km(parameters.gridsize); %deg2km(parameters.gridsize)*2;
parameters.raydensetol_azi=deg2km(parameters.gridsize_azi)*0.25; %deg2km(parameters.gridsize)*2;
parameters.fiterrtol = 2;   % error allowed in the wavelet fitting
parameters.dterrtol = 4;    % largest variance of the inversion error allowed
parameters.maxerrweight = 5; % Maximum error weight
parameters.polyfit_dt_err = 2; % (s) dt error greater than this, weighted 0