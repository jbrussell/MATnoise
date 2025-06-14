% Following Wapenaar et al. (2005) time-reversal approach.
%
% This simulates noise sources using true station geometry if 
% is_true_station_geometry=1. Otherwise, it uses a user defined 
% station geometry.
%
% Synthetic seismograms are saved as SAC files with the same structure
% as output from fetch_NOISE. These can then be used as input for
% the MATnoise ccf calculations.
%
% jbrussell - 6/2025

clear; close all;
rng default % for reproducibility

setup_parameters_synth;

% Path to output SAC files
path2out = ['./SAC/'];

is_true_station_geometry = 1; % Use real station geometry?

if is_true_station_geometry
    % Read real station geometry
    stalist = parameters.stalist;
    Y_stas = parameters.stalat;
    X_stas = parameters.stalon;
    z = parameters.staz;
    
    % Origin in lat lon
    olat = mean(Y_stas);
    olon = mean(X_stas);
else
    % Define custom station geometry
    % Origin in lat lon
    olat = 35;
    olon = -105;
    deg = 0.5; % station spacing

    % Build station geometry
    [Y_stas, X_stas] = meshgrid(olat+[-1:deg:1], olon+[-1:deg:1]);
    X_stas = X_stas(:);
    Y_stas = Y_stas(:);
    z = zeros(X_stas);
    % Perturb station locations so not a perfect grid (perfect grids cause artifacts in the beam pattern)
    d_deg = km2deg(25) * 0;
    X_stas = X_stas + d_deg*2*rand(length(X_stas),1)-d_deg;
    Y_stas = Y_stas + d_deg*2*rand(length(Y_stas),1)-d_deg;
    stalist = {};
    for ii = 1:length(Y_stas)
        stalist{ii} = ['STA',num2str(ii,'%.2d')];
    end
end


% Source timing
dt = 1; % [sec] sample rate
Ndays = 5; % Number of days to stack data
t = [0:dt:60*60*24-dt]; % [sec] time axis
tstart = datetime(2000,1,1); % initialize a datetime to start the experiment (it can be anything)

% Location of "source" ring (S)
N_sources = 1000; % Number of sources
N_excite = 1; % number of times to randomly excite each Source per hour
dr_S = km2deg(200); % [deg] width of annulus, S
r_S = km2deg(2000) + (dr_S*2*rand(N_sources,1)' - dr_S); %1000; % [km] radius of source ring
% dtheta = 0.1;
theta_S = 360*rand(N_sources,1)'; %[0:dtheta:360-dtheta];
% r_S(theta_S>90 & theta_S<180) = []; % remove sector of sources
% theta_S(theta_S>90 & theta_S<180) = []; % remove sector of sources
% amp_S = ones(size(theta_S)); % amplitude of sources
amp_S = ones(size(theta_S)) .* (cosd(theta_S-180+90+45)*5+6); % amplitude of sources

[y_S, x_S] = reckon(olat,olon,r_S,theta_S);

% Medium properties
vel = 3.5; % [km/s] velocity of medium

% Define source type
source_type = 'ricker'; % 'ricker' or 'microseism'
% RICKER
f_cent = 1/15; %1/8; % 1/100; % [1/s] dominant frequency of Ricker wavelet
% MICROSEISM
fmin = 1/10; % minimum frequency to sum over
fmax = 1/3; % maximum frequency to sum over

% =======================================================================

ccfpath = parameters.ccfpath;
nsta = parameters.nsta;

% Build frequency and time axes
Nt = length(t);
Fs = 1./dt;
f = Fs*(0:(Nt/2))/Nt;
time = ([0:Nt-1]-floor(Nt/2))*dt;  % build lagtime vector for plotting
time = [time(time<0), time(time>=0)];


figure(98); clf;
box on; hold on;
load coastlines
plot(coastlon,coastlat,'-b');
plot(X_stas,Y_stas,'o','color',[0.5 0.5 0.5],'linewidth',2)
% plot(x_Asrc,y_Asrc,'og','linewidth',2,'MarkerFaceColor','g');
% text(x_Asrc+50,y_Asrc,'A','fontsize',13);
% plot(x_Brec,y_Brec,'ob','linewidth',2,'MarkerFaceColor','b');
% text(x_Brec+50,y_Brec,'B','fontsize',13);
scatter(x_S,y_S,amp_S*10,amp_S,'o','filled');
text(x_S(1),y_S(1)+5,'S','fontsize',13);
axis square; axis equal;
set(gca,'fontsize',15,'linewidth',1.5,'layer','top');
xlabel('Lon (X)');
ylabel('Lat (Y)');
xlim([min(x_S)-1 max(x_S)+1]);
ylim([min(y_S)-1 max(y_S)+1]);
cb = colorbar;
ylabel(cb,'Source Amplitude');
colormap(viridis);

%% Load station data
R_Si_A = {};
for ista1 = 1:length(X_stas)
    sta1 = stalist{ista1};
        
    disp([sta1]);
        
    % Location of virtual source (A)
    x_Asrc = X_stas(ista1); %100; % km
    y_Asrc = Y_stas(ista1); %-100; % km
        

    %% Generate wavefield at ring of sources S

    % Calculate distance from Si's to A and B
    R_Si_A{ista1} = distance(y_S,x_S,y_Asrc,x_Asrc,referenceEllipsoid('GRS80'))/1000; % [km] distance from Si to A
    
    path2sta = [path2out,'/',sta1];
    if ~exist(path2sta)
        mkdir(path2sta);
    end
end
   

%% Generate synthetic dataset

for iday = 1:Ndays
    disp(['Day ',num2str(iday)])
    
    % Initialize wavefield at A due to random noise at all Si
    for ista1 = 1:length(X_stas)
        Si_A{ista1} = zeros(size(t));
    end
    
    for isrc = 1:length(x_S)

        % Loop through and generate sources, Si
        for ii = 1:N_excite
            % Generate random wavelet start times
            tshift = max(t) .* rand(1,1);
            % Generate random phase
            freq = f(f>=fmin & f<=fmax);
            phi_rand = (2*pi)*rand(1,length(freq));
            for ista1 = 1:length(Si_A)
                if strcmp(source_type,'ricker')
                    % Ricker wavelet
                    Si_A{ista1} = Si_A{ista1} + amp_S(isrc) .* ricker_wavelet(t-tshift,R_Si_A{ista1}(isrc),vel,f_cent);
                elseif strcmp(source_type,'microseism')
                    Si_A{ista1} = Si_A{ista1} + amp_S(isrc) .* microseism_source(t,R_Si_A{ista1}(isrc),vel,freq,phi_rand);
                else
                    error('Source type must be ''ricker'' or ''microseism''');
                end
            end
        end
    end
    
    for ista1 = 1:length(Si_A)
        sta1 = stalist{ista1};
        
        % Taper waveform
        Si_A{ista1} = cos_taper(Si_A{ista1});

        t_startsac = tstart+(iday-1);

        % If original SAC file does not exist, read required header info from the mat file
        H.KNETWK = 'XX'; %matin.traces(icomp).network;
        H.KSTNM = sta1; %matin.traces(icomp).station;
        H.KCMPNM = 'BHZ';%matin.traces(icomp).channel;
        H.KHOLE = 0; %matin.traces(icomp).location;
        H.STLA = Y_stas(ista1); %matin.traces(icomp).latitude;
        H.STLO = X_stas(ista1); %matin.traces(icomp).longitude;
        H.STEL = 0; %matin.traces(icomp).elevation;
        H.T0 = datenum(t_startsac); %matin.traces(icomp).startTime;

        % header origin time values
        tv = datevec(t_startsac);
        H.NZYEAR = tv(1);
        H.NZJDAY = datenum(tv(1:3)) - datenum(tv(1),1,1) + 1;
        H.NZHOUR = tv(4);
        H.NZMIN = tv(5);
        H.NZSEC = floor(tv(6));
        H.NZMSEC = (tv(6) - H.NZSEC)*1e3;

        H.DELTA = dt; %corrected.params.dt;
        H.NPTS = length(t); %length(corrected.params.taxis);

        data = Si_A{ista1}; %corrseis(corr_idx).timeseries;
        
%         {datadirectory}/{station}/{station}.{yyyy}.{jday}.{hh}.{mm}.{SS}.{COMP}.sac
        sac_path = [path2out,'/',sta1,'/',sta1,'.',num2str(H.NZYEAR),'.',num2str(H.NZJDAY,'%03d'),'.',num2str(H.NZHOUR,'%02d'),'.',num2str(H.NZMIN,'%02d'),'.',num2str(H.NZSEC,'%02d'),'.',H.KCMPNM,'.sac'];
%         sac_path = fullfile(sprintf('%s/%s.%s.%s.%s.sac',opath,evid, sacin.HEADER.KNETWK, sacin.HEADER.KSTNM, sacin.HEADER.KCMPNM));
        mksac(sac_path,data,H.T0,H); 
    end

end

%% Save source info

latS = y_S;
lonS = x_S;

save([ccfpath,'/sources.mat'],'N_sources','N_excite','dr_S','r_S','theta_S','amp_S','latS','lonS','x_S','y_S','vel','f_cent','Ndays','X_stas','Y_stas');

