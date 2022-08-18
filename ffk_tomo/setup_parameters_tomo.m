% Setup_parameters for raytomo
%
% https://github.com/jbrussell
addpath('../');
addpath('../functions');
addpath('./tomo_functions');
% parameters.qpath = './MINEOS_qfile/';

parameters.path2XSPsynth = []; %'./Xsp_synth/';
parameters.path2XSPsynth_SEM2D = []; %'./Xsp_synth_SEM2D/';

parameters.station_list = '../stations_all.txt';

% Lat Lon
parameters.lalim = [-9 -3] ;
parameters.lolim = [-136 -130];
parameters.gridsize = 0.1; % degrees
parameters.gridsize_azi = 0.1; %0.25; %3; %1.5; % gridsize for 2D azimuthal anisotropy (degrees)
parameters.r = 0.03; %0.01; % controls color bar [avgv(1-r) avgv(1+r)]

% RAYLEIGH 1st OVERTONE
parameters.comp = {'ZZ'}; % component
parameters.xspdir = 'ZZ_1S_LRT_N20'; % directory name
parameters.windir = 'window3hr_Zcorr_tiltcomp'; % window type
parameters.N_wl = 1; % number of wavelengths to consider during bessel fit measurements
parameters.frange = [1/25 1/3]; % [Hz]
parameters.per_ind = [2:2:20]; % index of periods to consider
% QC parameters
parameters.snr_tol = 3; % minimum signal-to-noise
parameters.r_tol_min = 0; % [km] minimum station separation
parameters.r_tol_max = 600; % [km] maximum station separation
parameters.err_tol = inf; % maximum misfit of bessel fit between observed and synthetic
parameters.is_rtolmin_wavelength = 1; parameters.wl_fac = 1.0; % determine distance tolerance by wavelength?
parameters.dep_tol = [0 0]; % [sta1 sta2] OBS depth tolerance
parameters.is_raydensity_thresh = 0; % Apply raydensity threshold to wipe out poorly constrained grid cells?
parameters.min_dep = 9999; %-3500 for min station depth to use

% Kernel Parameters
parameters.kernel_path = ['./SEM2D_FFK_save/',parameters.windir,'/']; % path to saved kernels
parameters.kmode = -1; % <0 : bandwidth; =0 instantaneous; >0 Fresnel zones
parameters.nfreqs = 30; % number of frequencies to include in bandwidth average
parameters.bw_frac = 0.25; % fraction of frequency to include in bandwidht
parameters.tphase_in = []; % Path to SEM traveltime surfaces. If left blank, will use analytical kernels by default.

% Smoothing parameters
parameters.damp0 = 1e-4; % Norm damping of isotropic phase velocity
parameters.smweight0 = 10; % isotropic second derivative smoothing
parameters.smweight0_azi = 1e3; %1000; % anisotropic second derivative smoothing
parameters.flweight0_azi = 0; %1000; % anisotropic first derivative flatness
parameters.damp0_azi = 0; % anisotropic norm damping (damp to zero)
parameters.is_wlsmooth = 0; %1; % weight smoothing by wavelength 0:no, 1:yes

% parameters for the tomography (QC)
parameters.raydensetol=deg2km(parameters.gridsize)*0.25; %deg2km(parameters.gridsize); %deg2km(parameters.gridsize)*2;
parameters.raydensetol_azi=deg2km(parameters.gridsize_azi)*0.25; %deg2km(parameters.gridsize)*2;
parameters.fiterrtol = 2;   % error allowed in the wavelet fitting
parameters.dterrtol = 4;    % largest variance of the inversion error allowed
parameters.maxerrweight = 5; % Maximum error weight
parameters.polyfit_dt_err = 2; % (s) dt error greater than this, weighted 0
parameters.stderrfac = 2; % remove measurements with err(i) > std(err)*stderrfac
