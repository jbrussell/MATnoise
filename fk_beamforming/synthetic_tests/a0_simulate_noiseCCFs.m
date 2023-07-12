% Following Wapenaar et al. (2005) time-reversal approach.
%
% This simulates noise sources using true station geometry if 
% is_true_station_geometry=1. Otherwise, it uses a user defined 
% station geometry.
%
% jbrussell - 7/2023

clear; close all;
rng default % for reproducibility

setup_parameters_synth;

is_overwrite = 0; %

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
Ndays = 5; %5; % Number of days to stack data
t = [0:dt:60*60]; % [sec] time axis

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
f_cent = 1/8; % 1/100; % [1/s] dominant frequency of Ricker wavelet
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

%% Generate station file

fid = fopen('stalist_synth.txt','w');
for ista = 1:length(X_stas)
    sta = stalist{ista};
    lat = Y_stas(ista);
    lon = X_stas(ista);
    fprintf(fid,'%10s %10f %10f %10f\n',sta,lat,lon,z(ista));
end
fclose(fid);

%% Generate synthetic dataset

for ista1 = 1:length(X_stas)
    sta1 = stalist{ista1};
    for ista2 = ista1:length(X_stas)
        
        sta2 = stalist{ista2};
        
        % Output path
        data_path = [ccfpath,'/window3hr/fullStack/ccfZZ/',sta1,'/'];
        if ~exist(data_path)
            mkdir(data_path);
        end
        if ~is_overwrite && exist([data_path,sta1,'_',sta2,'_f.mat'])
            disp([sta1,'-',sta2,' exists... skipping']);
            continue
        end
        
        disp([sta1,'-',sta2]);
        
        % Location of virtual source (A)
        x_Asrc = X_stas(ista1); %100; % km
        y_Asrc = Y_stas(ista1); %-100; % km

        % Location of receiver (B)
        x_Brec = X_stas(ista2); %-100; % km
        y_Brec = Y_stas(ista2); %100; % km
        
        figure(98); clf;
        box on; hold on;
        plot(coastlon,coastlat,'-b');
        plot(X_stas,Y_stas,'o','color',[0.5 0.5 0.5],'linewidth',2)
        plot(x_Asrc,y_Asrc,'og','linewidth',2,'MarkerFaceColor','g');
        text(x_Asrc,y_Asrc+2,'A','fontsize',13);
        plot(x_Brec,y_Brec,'ob','linewidth',2,'MarkerFaceColor','b');
        text(x_Brec,y_Brec+2,'B','fontsize',13);
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

        %% Generate wavefield at ring of sources S

        % Expected time between A and B
        R_A_B = distance(y_Asrc,x_Asrc,y_Brec,x_Brec,referenceEllipsoid('GRS80'))/1000; % [km] distance from A to B
        t_A_B = R_A_B ./ vel;

        % Theoretical prediction of ccf from Goran's class notes
        ccf_pre = real(1 ./ sqrt(t_A_B.^2 - time.^2));
        A = 1;
        J0_pre = besselj(0,2*pi*f*t_A_B) * A;

        % Theoretical prediction of spurious source arrival time
        [~,azi_A_B] = distance(y_Asrc,x_Asrc,y_Brec,x_Brec,referenceEllipsoid('GRS80')); % azi from A to B
        [~,azi_A_S] = distance(y_Asrc,x_Asrc,y_S(end),x_S(end),referenceEllipsoid('GRS80')); % azi from S to A
        azi_S_A = azi_A_S - 180;
        t_A_B_spurious = t_A_B .* cosd(azi_A_B-azi_S_A);
        
        % Calculate distance from Si's to A and B
        R_Si_A = distance(y_S,x_S,y_Asrc,x_Asrc,referenceEllipsoid('GRS80'))/1000; % [km] distance from Si to A
        R_Si_B = distance(y_S,x_S,y_Brec,x_Brec,referenceEllipsoid('GRS80'))/1000; % [km] distance from Si to B

        cohsum = zeros(size(t));
        ccf_auto = zeros(size(t));
        ihrs_total = 0;
        ccf_all = {};
        cohsum_all = {};
        ccf_misfit = [];
        for iday = 1:Ndays
            disp(['Day ',num2str(iday)])
            for ihr = 1 : (24/(max(t)/60/60))
%                 disp(ihr)
                % Generate wavefields at A and B due to random noise at all Si
                Si_A = zeros(size(t));
                Si_B = zeros(size(t));
                for isrc = 1:length(x_S)

                    % Loop through and generate sources, Si
                    for ii = 1:N_excite
                        if strcmp(source_type,'ricker')
                            % Ricker wavelet
                            % Generate random wavelet start times
                            tshift = max(t) .* rand(1,1);
                            Si_A = Si_A + amp_S(isrc) .* ricker_wavelet(t-tshift,R_Si_A(isrc),vel,f_cent);
                            Si_B = Si_B + amp_S(isrc) .* ricker_wavelet(t-tshift,R_Si_B(isrc),vel,f_cent);
                        elseif strcmp(source_type,'microseism')
                            freq = f(f>=fmin & f<=fmax);
                            phi_rand = (2*pi)*rand(1,length(freq));
                            Si_A = Si_A + amp_S(isrc) .* microseism_source(t,R_Si_A(isrc),vel,freq,phi_rand);
                            Si_B = Si_B + amp_S(isrc) .* microseism_source(t,R_Si_B(isrc),vel,freq,phi_rand);
                        else
                            error('Source type must be ''ricker'' or ''microseism''');
                        end
                    end
                end

                % Taper waveform
                Si_A = cos_taper(Si_A);
                Si_B = cos_taper(Si_B);

                % Calculate power spectra
                fftA = fft(Si_A);
                % fftA = spectrumwhiten_smooth(fftA,0.001);
        %         P_A = abs(fftA/Nt);
        %         P_A = P_A(1:Nt/2+1);
        %         P_A(2:end-1) = 2*P_A(2:end-1);
                fftB = fft(Si_B);
                % fftB = spectrumwhiten_smooth(fftB,0.001);
        %         P_B = abs(fftB/Nt);
        %         P_B = P_B(1:Nt/2+1);
        %         P_B(2:end-1) = 2*P_B(2:end-1);

                % Calculate cross-correlation (A --> B Causal; B --> A Acausal)
                cohsum = cohsum + fftB.*conj(fftA) ./ abs(fftB) ./ abs(fftA);
                ccf = real(ifft(cohsum,Nt)); % inverse FFT to get time domain
                ccf = fftshift(ccf); % rearrange values as [-lag lag]
                ccf = detrend(ccf);
                ccf = cos_taper(ccf);

                % Auto-correlation
                ccf_auto = ccf_auto + xcorr(Si_A,Si_A,floor(length(t)/2));

                ihrs_total = ihrs_total + 1;
                ccf_all{ihrs_total} = ccf;
                cohsum_all{ihrs_total} = cohsum;
                ccf_misfit(ihrs_total) = sum((ccf_all{ihrs_total} ./ max(ccf_all{ihrs_total}) - ccf_pre ./ max(ccf_pre)).^2) / length(ccf_pre);

                if 1 && mod(ihrs_total,Ndays*24)==0
                    figure(99); clf;
                    set(gcf,'position',[2         206        1730         814]);

                    subplot(2,2,1); box on; hold on;
                    plot(coastlon,coastlat,'-b');
                    plot(X_stas,Y_stas,'o','color',[0.5 0.5 0.5],'linewidth',2)
                    plot(x_Asrc,y_Asrc,'og','linewidth',2,'MarkerFaceColor','g');
                    text(x_Asrc,y_Asrc+2,'A','fontsize',13);
                    plot(x_Brec,y_Brec,'ob','linewidth',2,'MarkerFaceColor','b');
                    text(x_Brec,y_Brec+2,'B','fontsize',13);
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

                    subplot(2,2,3); box on; hold on;
                    plot([1:ihrs_total],ccf_misfit,'-b','linewidth',2);
                    plot(ihrs_total,ccf_misfit(ihrs_total),'ok','linewidth',1,'markersize',12,'MarkerFaceColor','b');
                    set(gca,'fontsize',15,'linewidth',1.5);
                    xlabel('Data length (hr)');
                    ylabel('Misfit');

        %             % Filter CCF before plotting
        %             coperiod = [5 20];
        %             costap_wid = 0.2; % 0 => box filter; 1 => Hann window
        %             [ ccf_filtered ] = tukey_filt( fft(fftshift(ccf_all{ihrs_total})),coperiod,dt,costap_wid );
        %             ccf_all{ihrs_total} = fftshift(real(ifft(ccf_filtered)));

                    subplot(2,2,2); box on; hold on;
        %             ccf_pre_conv = conv(ccf_pre,ccf_auto,'same');
                    plot(time,ccf_all{ihrs_total} ./ max(ccf_all{ihrs_total}),'-k','linewidth',2);
                    plot(time,ccf_pre ./ max(ccf_pre),'-r','linewidth',1);
        %             plot(time,ccf_pre_conv ./ max(ccf_pre_conv),'--g','linewidth',1);
                    plot(t_A_B_spurious*[1 1],[0 1],'--r','linewidth',1);
                    set(gca,'fontsize',15,'linewidth',1.5)
                    xlim([-200 200]);
                    title([num2str(ihrs_total),' hrs']);
                    xlabel('Lag Time (s)');
                    ylabel('CCF_{A,B}');

                    P_ccf = abs(cohsum_all{ihrs_total}/Nt);
                    P_ccf = P_ccf(1:Nt/2+1);
                    P_ccf(2:end-1) = 2*P_ccf(2:end-1);

        %             subplot(2,1,2); box on; hold on;
        %             plot(1./f,P_ccf,'-b');
        %             plot(1./f,smooth(P_ccf,20),'-k','linewidth',2);
        %             plot(1./f_cent*[1 1],[min(P_ccf) max(P_ccf)],'--r','linewidth',1.5);
        %             % xlim([0 500]);
        %             set(gca,'fontsize',15,'linewidth',1.5,'xscale','log','yscale','log')
        %             xlabel('Period (s)');

                    subplot(2,2,4); box on; hold on;
                    J0 = real(cohsum_all{ihrs_total}(1:Nt/2+1));
                    plot(f,J0 ./ rms(J0),'-b');
                    plot(f,smooth(J0,10) ./ rms(J0),'-k','linewidth',2);
                    plot(f,J0_pre ./ rms(J0_pre),'-r','linewidth',1);
                    plot(f_cent*[1 1],[0 1],'--r','linewidth',1.5);
                    % xlim([0 500]);
                    set(gca,'fontsize',15,'linewidth',1.5,'xscale','linear','yscale','linear')
                    xlabel('Freq (Hz)');
                    ylabel('J_0');
                    xlim([0 max(f)]);

                    drawnow;
                end
            end
        end
        
        % Save data
        lon1 = x_Asrc; 
        lat1 = y_Asrc; 
        lon2 = x_Brec; 
        lat2 = y_Brec; 
        coh_sum = cohsum_all{end};
        coh_num = ihrs_total;
        coh_sum_win = coh_sum;
        ccf = ccf_all{end};
        min_grv = [];
        max_grv = [];
        stapairsinfo.stanames = {sta1, sta2};
        stapairsinfo.lats = [lat1, lat2];
        stapairsinfo.lons = [lon1, lon2];
        stapairsinfo.dt = dt;
        stapairsinfo.r = R_A_B;
        save([data_path,sta1,'_',sta2,'_f.mat'],'coh_num','coh_sum','coh_sum_win','max_grv','min_grv','stapairsinfo','ccf','time');
        
        % Save reverse station pair
        if ista1~=ista2
            coh_sum = conj(cohsum_all{end});
            coh_num = ihrs_total;
            coh_sum_win = coh_sum;
            ccf = flip(ccf_all{end});
            min_grv = [];
            max_grv = [];
            stapairsinfo.stanames = {sta2, sta1};
            stapairsinfo.lats = [lat2, lat1];
            stapairsinfo.lons = [lon2, lon1];
            stapairsinfo.dt = dt;
            stapairsinfo.r = R_A_B;
            data_path = [ccfpath,'/window3hr/fullStack/ccfZZ/',sta2,'/'];
            if ~exist(data_path)
                mkdir(data_path);
            end
            save([data_path,sta2,'_',sta1,'_f.mat'],'coh_num','coh_sum','coh_sum_win','max_grv','min_grv','stapairsinfo','ccf','time'); 
        end
        
    end

end

%% Save source info

latS = y_S;
lonS = x_S;

save([ccfpath,'/sources.mat'],'N_sources','N_excite','dr_S','r_S','theta_S','amp_S','latS','lonS','x_S','y_S','vel','f_cent','Ndays','X_stas','Y_stas');

