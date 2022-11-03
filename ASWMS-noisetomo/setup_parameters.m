% Script to setup parameters used for the whole project

% prepend functions directory to MATLAB path
fullMAINpath = mfilename('fullpath');
functionspath = [fullMAINpath(1:regexp(fullMAINpath,mfilename)-1),'functions'];
plottingpath = [fullMAINpath(1:regexp(fullMAINpath,mfilename)-1),'plotting_scripts'];
addpath(functionspath); addpath(plottingpath);

parameters.workingdir = './OUT/';

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
%%%% Global settings
parameters.proj_name = 'XX'; % network code
parameters.component = 'BHZ';   % determined by filenames
parameters.lalim = [-9 -3] ;
parameters.lolim = [-136 -130];
parameters.gridsize = 0.25; % degrees?

% Parameters for reading in matnoise Bessel fits
% RAYLEIGH 1st OVERTONE
parameters.station_list = './stations_all.txt';
parameters.compnoise = {'ZZ'}; % component
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

if exist('periods.mat')
    load('periods.mat');
    parameters.periods = periods;
else
%     parameters.periods = round(logspace(log10(20),log10(150),15));
    parameters.periods = ones(1,length(parameters.per_ind));
end


%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% % parameters for data downloading (if using IRIS DMC)
% parameters.start_time = '2009-01-07 00:00:00';
% parameters.end_time = '2009-06-08 00:00:00'; % put '' for using 4 days before current date
% parameters.is_use_timestamp = 0;
% parameters.network = '_US-ALL';
% parameters.minMw = 6;
% parameters.maxdepth = 50;
% parameters.datalength = 7200;  % in second
% parameters.resample_delta = 1; % in second
% parameters.dbpath = './sacdata/';
% parameters.eventfile = 'eventlist';

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% Parameters for own data selection criteria (if using SAC)
parameters.dbpath = '/path/to/sac/data/directory/'; 
parameters.eventfile = 'evlist.txt'; % name of eventlist located in dbpath
parameters.minMw = 5.5;
parameters.maxdepth = 50;
parameters.snr_tol = 3;
parameters.resample_delta = 1; % in second 

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for the auto_win_select.m
parameters.largest_epidist_range = 3000; % [km] maximum array size to consider
parameters.cycle_before = 2; % number of cycles to buffer at start of window
parameters.cycle_after = 5; % number of cycles to buffer at end of window
parameters.min_dist_tol = deg2km(20); % min epicentral distance
parameters.max_dist_tol = deg2km(100); %deg2km(160); % min epicentral distance
parameters.min_groupv = 2; % minimum group velocity as a starting guess
parameters.max_groupv = 5; % maximum group velocity as a starting guess
parameters.cent_freq = mean(1./parameters.periods); %0.025;
parameters.min_sta_num = 5; %10 %JBR
parameters.is_removebadwindows = 1; % JBR - 1:remove bad windows (default), 0:keep all windows (This option may be desirable if lacking data. Will need to manually adjust windows using run_exam_winpara.m)

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for the cross-correlation measurement
% (gsdfmain.m)
parameters.minstadist = 5; % minimum station cross-correlation distance in km
parameters.maxstadist = 600; %600; %250 %200;   % maximum station cross-correlation distance in km
parameters.is_rm_resp = 0; % remove instrument response?
parameters.periods = sort(parameters.periods);  % make sure periods are ascending
parameters.refv = 4;   % to select the correct cycle (for cycle skipping correction)
parameters.refphv = ones(size(parameters.periods))*4;
parameters.min_width = 0.10; %0.06;  % to build up gaussian filters
parameters.max_width = parameters.min_width; %0.10;  
parameters.wintaperlength = 30;   % taper to build up the isolation filter
parameters.prefilt_frac = 0.25; % fraction of frequency beyond min and max to prefilter to
parameters.prefilter = [1/(1/min(parameters.periods)*(1+parameters.prefilt_frac)) 1/(1/max(parameters.periods)*(1-parameters.prefilt_frac))]; %[15,75]; %[15,160]; %[10,200];
parameters.xcor_win_halflength = 1.5*max(parameters.periods); %200; %300; %200 %150;  % window length for the cross-correlation
parameters.xcor_win_iter = zeros(size(parameters.periods)); % re-apply the xcor window due to measured group delay, should be same length as periods, not used anymore
parameters.Nfit = 2; %4; % 2 % number of cycles fit on either side of maximum for gaussian wavelet
parameters.Ncircle = 5; % number of cycles searched for cycle skipping
parameters.cohere_tol = 0.65; % minimum coherenecy between two stations
parameters.tp_tol = 10;  % seconds away from averaged phase velocity 
parameters.is_winx2 = 1; % jbr - 0=window single seismogram (default); 1=window both seismograms (may help reduce overtone interference)

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for the tomography
% (eikonal_eq.m helmholtz_eq.m)
% parameters.smweight_array = 3*[0.4 0.3 0.2 0.2 0.2 0.5 0.5 0.5 1 1 1 2 3 3 3];  % 2nd derivative smoothing weight for the deltaSx and delta Sy
parameters.grd_per_wl = parameters.refphv.*parameters.periods./deg2km(parameters.gridsize); % JBR - grid cells per wavelength
parameters.smweight_array = 100*0.1*parameters.grd_per_wl; % 2nd derivative smoothing weight for the deltaSx and delta Sy
parameters.flweight_array = 0*parameters.grd_per_wl; % JBR - 1st derivative smoothing
parameters.is_offgc_smoothing = 0; % 1st derivative smoothing along propagation direction rather than great circle. Requires an initial run of a6_a0_eikonal_eq_GetPropAzi.m to get propagation azimuth
parameters.raydensetol=deg2km(parameters.gridsize)*2;
parameters.Tdumpweight = 0;  % damping the ray to the great-circle path
parameters.Rdumpweight = 0;  % damping the region to have the same phase velocity
parameters.fiterrtol = 3;   % error allowed in the wavelet fitting
parameters.isRsmooth = 1;  % smoothing due to Sx and Sy or Sr and S_theta; 1=Sr,S_theta
parameters.dterrtol = 2;    % largest variance of the inversion error allowed
parameters.inverse_err_tol = 2;  % count be number of standard devition
parameters.min_amp_tol = 0.1;  % station with amplitude smaller than this ratio of average amplitude will not be used.
parameters.amp_var_tol = 0.4; % station with amplitude +/- this fraction relative to mean value of nearby stations will not be used (smaller value, more restrictive)
parameters.alpha_range = [1 1];
parameters.alpha_search_grid = 0.1;
parameters.Nwl_mask = 2; % number of wavelengths under which to mask

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameter for stacking 
% (stack_phv.m stack_helm.m)
parameters.min_csgoodratio= 1*ones(1,length(parameters.periods));%[3 3 3 3 5 10 15 15 15 15 20]; %[3 3 3 3 5 10 15 15]; % minimum radio between good and bad measurements for a good event
parameters.min_phv_tol = 3; % minimum allowed phase velocity
parameters.max_phv_tol = 5; % maximum allowed phase velocity
parameters.is_raydense_weight = 0; %1; % manual says suggested turned off for large azimuthal anisotropy
parameters.min_event_num = 3; %10;
parameters.err_std_tol = 2;
parameters.issmoothmap = 1;
parameters.smooth_wavelength = 0.25;
parameters.event_bias_tol = 3; %2;


%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for azimuthal anisotropy inversion
% beta version
parameters.gridsize_azi = parameters.gridsize; % grid size for azimuthal anisotropy in a6_b_eikonal_2DanisoRT.m
parameters.smsize = 3; %1;  % averaging nearby grid number
parameters.azi_bin_deg_ani= 20; % [deg] size of azimuthal bins
parameters.off_azi_tol = 30; % differ from great circle path in degrees
parameters.is_one_phi = 0; %1;
parameters.is_offgc_propagation = 0; % % Account for off-great-circle propagation using eikonal tomography maps? Otherwise will assume great-circle propagation.

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% Attenuation
% parameters for amplitude/gradient field estimates and receiver terms
parameters.is_receiver_terms = 0; % Correct amplitudes using receiver terms calculated from a8_0_receiver_terms?
parameters.is_eikonal_ampgrad = 1; % 1: use eikonal tomography values for amplitude gradient; 0: use amplitude field estimates
parameters.is_eikonal_ampgrad_norm = 1; % 1: use grad(A)/A values from direct inversion, forgoing the need to independently determine the amplitude field.
parameters.is_eikonal_phasegrad = 1; % 1: use eikonal tomography values for phase gradient; 0: use travel-time field estimates

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for solving for receiver terms
% (a8_0_receiver_terms)
parameters.max_sta_dist = 150; % [km] maximum separation allowed between station pairs
parameters.is_azibin = 1; % bin data by propagation azimuth?
parameters.deg_bins = 15; % [deg] size of azimuthal bins in degrees
parameters.avg_type = 'median'; % 'median'; 'mean'
parameters.min_Mw_rec = 5.5; %6.0;

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for solving for estimating attenuation coefficient alpha and
% amplification beta
% (b1_estimate_alpha_beta)
parameters.min_Mw_alpha = 5.5; % minimum magnitude
parameters.min_Ngrcells = 20; % minimum numbe of grid cells required in order to use event
parameters.azi_bin_deg = 20; % [deg] size of azimuthal bins
parameters.min_nbin = 10; % minimum number of measurements in order to include bin
parameters.N_min_evts = 10; % minimum number of events contributing to grid cell required in order to be considered
parameters.smsize_alpha = 3; % number of nearby gridcells to gather data from
parameters.smweight_beta = 0.3; % Second derivative smoothing weight for beta map
parameters.smooth_alpha_Nwl = 1.5; % [wavelengths] smoothing radius of 2d alpha map
parameters.azi_anom_maxthresh = 10; % [degrees] Remove grid cells with propagation anomaly larger than this value
parameters.isbin_2d = 1; % use azimuthal binning for 2-D maps?

%%
if length(parameters.periods)~=length(parameters.smweight_array) || length(parameters.periods)~=length(parameters.min_csgoodratio)
    error('Length of periods doesn''t match smweight_array and/or min_csgoodratio');
end

system(['cp ./setup_parameters.m ',parameters.workingdir]);
