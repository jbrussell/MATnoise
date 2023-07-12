% Run beamforming roughly following Gal et al. (2014) doi: 10.1093/gji/ggu183 
% "Improved implementation of the fk and Capon methods for array analysis 
% of seismic noise"
%
% Modeled code after: https://geophydog.cool/post/ncf_cross_spectral_beamforming/
%
% jbrussell - 7/2023

clear; close all;

setup_parameters_synth;

%======================= PARAMETERS =======================%
comp = {'ZZ'}; %'RR'; 'ZZ'; 'TT'
windir = 'window3hr';

% Data weights and SNR threshold
is_azimuthal_weight = 1; % apply azimuthal weights to downweight redundant azimuths?
dazi = 10; % [deg] bin width for azimuthal homogenization
is_dist_weight = 1; % down weight shorter station separations
snr_thresh = 3; % minimum SNR to consider
min_grv_snr = 1.2; % minimum group velocity for SNR window
max_grv_snr = inf; % maximum group velocity for SNR window

is_all_freq = 0; % Use all available frequencies between per_min and per_max? (takes longer)

% Time window to consider
t_min = -200;
t_max = 200;

% Periods to average over
per_min = 5; % [sec] minimum period
per_max = 10; % [sec] maximum period
Npers = 30; % number of periods to consider
f_vec = linspace(1/per_max,1/per_min,Npers);
per_vec = 1./f_vec;

% Slowness values to search over
s_min = 1/5; % s/km
s_max = 1/2.5; % s/km
Nslow = 100; % s/km
s_vec = linspace(s_min,s_max,Nslow);

% Back-azimuth values to search over
Nbaz = 360;
baz_vec = linspace(0,360,Nbaz);

[s_mat, baz_mat] = meshgrid(s_vec,baz_vec);

%%%%%

ccfpath = parameters.ccfpath;

%% If synthetic, plot source geometry
synth_sources = [ccfpath,'/sources.mat'];
if exist(synth_sources)
    load(synth_sources);
    
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
    drawnow
end

%% Determine data weights based on azimuth (downweight common azimuths)
% input path
ccf_path = [ccfpath,windir,'/fullStack/ccf',comp{1},'/'];

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
        [snr, signal_ind] = calc_SNR(ccf_filtered,min_grv_snr,max_grv_snr,r1,isfigure_snr);
        
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

%% Do beamforming

%%% --- Loop through station 1 --- %%%
Pf = zeros(Nbaz,Nslow,Npers);
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
        if is_all_freq
            freq = f(f>=1/per_max & f<=1/per_min);
            per_vec = 1./freq;
            Npers = length(per_vec);
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

%% Plot beam

P_abs = abs(P);
P_abs = P_abs / max(P_abs(:));
P_abs = 10*log10(P_abs);
x = s_mat.*sind(baz_mat);
y = s_mat.*cosd(baz_mat);

figure(1); clf;
set(gcf,'position',[616   587   504   431]);
[h,c] = polarPcolor(s_vec,baz_vec,P_abs,'Nspokes',9,'fontsize',13);
colormap(viridis);
c.LineWidth = 1.5;
ylabel(c,'Relative Power (dB)','fontsize',15);
set(gca,'fontsize',15,'linewidth',1.5)
caxis([prctile(P_abs(:),80) 0]);
titl = title([num2str(per_min),'-',num2str(per_max),'s']);
titl.Position(2) = titl.Position(2) + 0.25;
if exist(synth_sources)
    hp = polar([0:360]*pi/180,1/vel/(s_max-s_min)*ones(size(0:360)),'--r');
    hp.LineWidth = 1.5;
end


