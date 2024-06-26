% Run beamforming roughly following Gal et al. (2014) doi: 10.1093/gji/ggu183 
% "Improved implementation of the fk and Capon methods for array analysis 
% of seismic noise"
%
% Modeled code after: https://geophydog.cool/post/ncf_cross_spectral_beamforming/
%
% jbrussell - 7/2023
%
% Update 6/20: Add functionality to apply a group velocity window in the time
% domain. This allows for separation of different mode branches at short
% periods in oceanic setting.
%
% Update 6/20: Assume slowness is known within a small range and produce plots
% of power versus back azimuth. In practice, we do this by calculating
% all slowness values in setup_parameters_beamforming then stacking power only 
% over the range of assumed slowness values.

clear; close all;

setup_parameters_beamforming;

is_save_mat = 1; % save output mat file?

%======================= PARAMETERS =======================%
is_azimuthal_weight = 0; % apply azimuthal weights to downweight redundant azimuths?
dazi = 10; % [deg] bin width for azimuthal homogenization
is_dist_weight = 0; % down weight shorter station separations
snr_thresh = 3; % minimum SNR to consider
min_grv_snr = 1.2; % minimum group velocity for SNR window
max_grv_snr = inf; % maximum group velocity for SNR window

is_grv_window = 1; % apply a group velocity window for analysis? 1 or 0
min_grv_win = 1.6; % [km/s] minimum group velocity for analysis window
max_grv_win = inf; % [km/s] maximum group velocity for analysis window

s_min_assumed = 1/4.0; % [s/km] minimum of assumed slowness range
s_max_assumed = 1/3.9; % [s/km] maximum of assumed slowness range

is_all_freq = 0; % Use all available frequencies between per_min and per_max? (takes longer)

%%%%
comp = parameters.comp;
windir = parameters.windir;

% Time window to consider
t_min = parameters.t_min;
t_max = parameters.t_max;

% Periods to average over
per_min = parameters.per_min;
per_max = parameters.per_max;
Npers = parameters.Npers;
f_vec = linspace(1/per_max,1/per_min,Npers);
per_vec = 1./f_vec;

% Slowness values to search over
s_min = parameters.s_min;
s_max = parameters.s_max;
Nslow = parameters.Nslow;
s_vec = linspace(s_min,s_max,Nslow);

% Back-azimuth values to search over
Nbaz = parameters.Nbaz;
baz_vec = linspace(0,360,Nbaz);

%% Plot station geometry
figpath = parameters.figpath;
if ~exist(figpath)
    mkdir(figpath)
end

figure(98); clf;
set(gcf,'color','w');
box on; hold on;
load coastlines
plot(coastlon,coastlat,'-b');
plot(stalon,stalat,'ok','markerfacecolor','r','markersize',8,'linewidth',1.5)
axis square; axis equal;
set(gca,'fontsize',15,'linewidth',1.5,'layer','top');
xlabel('Lon (X)');
ylabel('Lat (Y)');
xlim([min(stalon)-1 max(stalon)+1]);
ylim([min(stalat)-1 max(stalat)+1]);
drawnow
    
save2pdf([figpath,'array_geom.pdf'],98,250);

%% Determine data weights based on azimuth (downweight common azimuths)
% input path
ccf_path = [parameters.ccfpath,windir,'/fullStack/ccf',comp{1},'/'];

stalist = parameters.stalist;
nsta=parameters.nsta; % number of target stations to calculate for

% Loop through all stations and gather info
stainfo = [];
ii = 0;
for ista1= 1:nsta    
    sta1=char(stalist(ista1,:));
    sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
    for ista2 = 1: nsta 
        sta2 = char(stalist(ista2,:));
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        
        
        %%% --- Load in the ccf --- %%%
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        
        if ~exist(filename,'file')
            disp(['not exist ',filename])
            continue;
        end
        disp([sta1,'-',sta2]);
        
        data1 = load(filename);
        r1 = distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2),referenceEllipsoid('GRS80'))/1000;
        [~,az] = distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2),referenceEllipsoid('GRS80'));
        
        ccf_fft = data1.coh_sum ./ data1.coh_num;
        ccf_fft(isnan(ccf_fft)) = 0;
        
        %----------- Frequency ==> Time domain -------------%
        N = length(ccf_fft);
        ccf = real(ifft(ccf_fft,N)); % inverse FFT to get time domain
        ccf = fftshift(ccf); % rearrange values as [-lag lag]
        ccf = detrend(ccf);
        ccf = cos_taper(ccf);
        %----------- FILTER DATA (FREQUENCY DOMAIN) -------------%
        costap_wid = 0.2;
        dt = data1.stapairsinfo.dt;
        [ ccf_filtered ] = tukey_filt( fft(fftshift(ccf)),[per_min per_max],dt,costap_wid );
        ccf = fftshift(real(ifft(ccf_filtered)));
        %----------- SIGNAL TO NOISE RATIO -------------%
        isfigure_snr = 0;
        [snr, signal_ind] = calc_SNR(ccf_filtered,min_grv_snr,max_grv_snr,r1,dt,isfigure_snr);
        
        if snr < snr_thresh
            continue
        end
        
        ii = ii + 1;
        stainfo.r(ii) = r1;
        stainfo.az(ii) = az;
        stainfo.sta1{ii} = sta1;
        stainfo.sta2{ii} = sta2;
        stainfo.snr(ii) = snr;
    end
end

stainfo.w = ones(size(stainfo.r));

% Down weight redundant azimuths
if is_azimuthal_weight
    az_bins = [0:10:360];
    az_weights = ones(size(az_bins));
    for ibin = 1:length(az_bins)-1
        I = find(stainfo.az>=az_bins(ibin) & stainfo.az<az_bins(ibin+1));
        az_weights(ibin) = 1/length(I);
        stainfo.w(I) = stainfo.w(I) * (1/length(I));
    end
end

% Down weight shorter station separations
if is_dist_weight
    dist_weights = (stainfo.r/max(stainfo.r));
    stainfo.w = stainfo.w .* (stainfo.r/max(stainfo.r));
%     stainfo.w(stainfo.r<200) = 0;
end

%% Load data and do beamforming

%%% --- Loop through station 1 --- %%%
Pf = zeros(Nbaz,Nslow,Npers);
itestflag = 0;
for ista1= 1:nsta
    
    sta1=char(stalist(ista1,:));
    sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
    
    %%% --- Loop through station 2 --- %%%
    for ista2 = 1: nsta % length(v_sta)
        sta2 = char(stalist(ista2,:));
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
                
        %%% --- Load in the ccf --- %%%
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        
        if ~exist(filename,'file')
            disp(['not exist ',filename])
            continue;
        end
        disp([sta1,'-',sta2]);
        
        data1 = load(filename);
        r1 = distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2),referenceEllipsoid('GRS80'))/1000;
        [~,az] = distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2),referenceEllipsoid('GRS80'));
        
        % Get weight for station pair
        Ipair = find(strcmp(stainfo.sta1,sta1) & strcmp(stainfo.sta2,sta2));
        w = stainfo.w(Ipair);
        if isempty(w)
            continue
        end
        
        % Displacement vector pointing from station 1 to station 2
        %          x             y
        r_v = r1*[sind(az); cosd(az)];
        
        % ccf in frequency domain
        ccf_fft = data1.coh_sum ./ data1.coh_num;
        
        % Build frequency and time axes
        Nt = length(ccf_fft);
        dt = data1.stapairsinfo.dt;
        Fs = 1./dt;
        f = Fs*(0:(Nt/2))/Nt; % positive frequencies
        time = ([0:Nt-1]-floor(Nt/2))*dt;  % build lagtime vector for plotting
        time = [time(time<0), time(time>=0)];
        
        % Do group velocity window
        if is_grv_window
            isfigure = 1;
            [ccf_fft_win] = window_ccf(ccf_fft, min_grv_win, max_grv_win, r1, dt, isfigure);
            ccf_fft = ccf_fft_win;
        end
        
        % Index desired time window
        ccf = real(ifft(ccf_fft,Nt)); % inverse FFT to get time domain
        ccf = fftshift(ccf); % rearrange values as [-lag lag]
        ccf = detrend(ccf);
        ccf = cos_taper(ccf);
        ccf(time<t_min | time>t_max) = 0;
        ccf(time>=t_min & time<=t_max) = cos_taper(ccf(time>=t_min & time<=t_max));
        % Shift back to frequency domain
        ccf_fft = fft(fftshift(ccf));
        
        % Index positive frequencies
        ccf_fft = 2*ccf_fft(1:Nt/2+1);
        
        % Use all available frequencies between per_min and per_max
        if is_all_freq && itestflag==0
            freq = f(f>=1/per_max & f<=1/per_min);
            per_vec = 1./freq;
            Npers = length(per_vec);
            Pf = zeros(Nbaz,Nslow,Npers);
            itestflag = 1;
        end
        for iper = 1:Npers
            per = per_vec(iper);
            omega = 2*pi ./ per;
            ccf_fft_i = interp1(f,ccf_fft,1/per);
            
            for islow = 1:Nslow
                slow = s_vec(islow);
                
                for ibaz = 1:Nbaz
                    baz = baz_vec(ibaz);
                    
                    % Slowness vector pointing from source N toward array
                    %                  x             y
                    s_v = slow*[sind(baz-180); cosd(baz-180)]; 
                    
                    % slowness direction dotted onto interstation path
                    sdotr = s_v'*r_v;
                    
                    % Estimate power
                    Pf(ibaz,islow,iper) = Pf(ibaz,islow,iper) + w * ccf_fft_i * exp(1i*omega * sdotr);

                end
            end
        end
%         figure(1); clf;
%         % scatter(x(:),y(:),10,P_abs(:));
%         % contourf(x,y,P_abs,'LineColor','none'); axis equal;
%         % polarscatter(baz_mat(:)*pi/180,s_mat(:),10,P_abs(:));
%         [h,c] = polarPcolor(s_vec,baz_vec,abs(P),'Nspokes',9);
%         colormap(viridis);
    end
end

% Sum over frequencies
P = sum(Pf,3);

% Sum over assumed slowness values
Islow = find(s_vec>=s_min_assumed & s_vec<=s_max_assumed);
P_s = sum(P(:,Islow),2);

%% Plot beam

P_abs = abs(P);
P_abs = P_abs / max(P_abs(:));
P_abs = 10*log10(P_abs);

figure(1); clf;
subplot(2,1,1);
set(gcf,'position',[616   130   562   888],'color','w');
[h,c] = polarPcolor(s_vec,baz_vec,P_abs,'Nspokes',9,'fontsize',13);
colormap(viridis);
c.LineWidth = 1.5;
ylabel(c,'Relative Power (dB)','fontsize',15);
set(gca,'fontsize',15,'linewidth',1.5)
caxis([prctile(P_abs(:),80) 0]);
titl = title([num2str(per_min),'-',num2str(per_max),'s']);
titl.Position(2) = titl.Position(2) + 0.1;
titl.Position(1) = titl.Position(1) + 0;

% Plot bounds of assumed slowness
hp = polar([0:360]*pi/180,s_min_assumed/(s_max-s_min)*ones(size(0:360)),'-r');
hp.LineWidth = 1.5;
hp = polar([0:360]*pi/180,s_max_assumed/(s_max-s_min)*ones(size(0:360)),'-r');
hp.LineWidth = 1.5;

% Plot power by azimuth
P_s_abs = abs(P_s);
P_s_abs = P_s_abs / max(P_s_abs(:));
P_s_abs = 10*log10(P_s_abs);
subplot(2,1,2);
plot(baz_vec,P_s_abs,'-r','linewidth',3);
xlabel('Back Azimuth (\circ)');
ylabel('Relative Power (dB)');
title(['Stacked between ',num2str(s_min_assumed),'-',num2str(s_max_assumed),' s/km'])
set(gca,'fontsize',15,'linewidth',1.5)
xlim([0 360]);

save2pdf([figpath,'fk_beamform_',comp{:},'_',num2str(per_min),'_',num2str(per_max),'s_grvwindow_',num2str(min_grv_win),'_',num2str(max_grv_win),'kms_phv_',num2str(1/s_max_assumed),'_',num2str(1/s_min_assumed),'kms.pdf'],1,250);

%% Save beam

matpath = './beam_out/';
if is_save_mat
    if ~exist(matpath)
        mkdir(matpath)
    end
    
    beam.P = P;
    beam.P_abs = P_abs;
    beam.s_vec = s_vec;
    beam.baz_vec = baz_vec;
    beam.parameters = parameters;

    save([matpath,'fk_beamform_',comp{:},'_',num2str(per_min),'_',num2str(per_max),'s_grvwindow_',num2str(min_grv_win),'_',num2str(max_grv_win),'kms_phv_',num2str(1/s_max_assumed),'_',num2str(1/s_min_assumed),'kms.mat'],'beam');
end


