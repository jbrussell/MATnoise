% Calculate Array Response Function (ARF). The ARF depends only on the 
% geometry of the seismic array and not on the data, but running it
% after the ccfs have been calculated allows us to examine the effect of
% SNR threshold and distance and azimuth weighting as would be done
% with the real dataset. 
%
% See Eq. 2.16 of "Seismic Ambient Noise" by Nakata, Gualtieri, & Fichtner
%
% jbrussell - 2/2024

clear; close all;

setup_parameters_beamforming;

%========== Define plane wave properties for Array Response Function ==========%

f_p = 1/8; % [Hz] frequency of plane wave
az_p = 0; %[deg] propagation azimuth. 0=North propagation (i.e., arriving from south)
slow_p = 0; %1/3.5; %[s/km] slowness of plane wave (zero gives vertically incident wave)
% Slowness vector pointing from source N toward array
%                  x           y
s_p = slow_p*[sind(az_p); cosd(az_p)]; 

%======================= PARAMETERS =======================%
is_azimuthal_weight = 1; % apply azimuthal weights to downweight redundant azimuths?
dazi = 10; % [deg] bin width for azimuthal homogenization
is_dist_weight = 1; % down weight shorter station separations
snr_thresh = 3; % minimum SNR to consider
min_grv_snr = 1.2; % minimum group velocity for SNR window
max_grv_snr = inf; % maximum group velocity for SNR window

is_all_freq = 0; % Use all available frequencies between per_min and per_max? (takes longer)

%%%%
comp = parameters.comp;
windir = parameters.windir;


% Periods to average over
per_min = 1/f_p; %5; % [sec] minimum period
per_max = 1/f_p; %10; % [sec] maximum period
Npers = 1;  %30; % number of periods to consider
f_vec = linspace(1/per_max,1/per_min,Npers);
per_vec = 1./f_vec;

% Slowness values to search over
s_min = 0; %1/5; % s/km
s_max = 1/2.5; % s/km
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
        
        for iper = 1:Npers
            per = per_vec(iper);
            omega = 2*pi ./ per;
            
            for islow = 1:Nslow
                slow = s_vec(islow);
                
                for ibaz = 1:Nbaz
                    baz = baz_vec(ibaz);
                    
                    % Slowness vector pointing from source N toward array
                    %                  x             y
                    s_v = slow*[sind(baz-180); cosd(baz-180)]; 
                    
                    % slowness direction dotted onto interstation path
                    sdotr = (s_v-s_p)'*r_v;
                    
                    % Estimate power
                    Pf(ibaz,islow,iper) = Pf(ibaz,islow,iper) + w * exp(1i*omega * sdotr);
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

figure(1); clf;
set(gcf,'position',[616   587   504   431],'color','w');
[h,c] = polarPcolor(s_vec,baz_vec,P_abs,'Nspokes',9,'fontsize',13);
colormap(viridis);
c.LineWidth = 1.5;
ylabel(c,'Relative Power (dB)','fontsize',15);
set(gca,'fontsize',15,'linewidth',1.5)
caxis([prctile(P_abs(:),80) 0]);
titl = title([num2str(per_min),'-',num2str(per_max),'s']);
titl.Position(2) = titl.Position(2) + 0.25;
hp = polar((az_p+180+90)*pi/180,slow_p/(s_max-s_min),'-or');
hp.LineWidth = 1.5;
hp.MarkerEdgeColor = 'w';
hp.MarkerFaceColor = 'r';
sg = sgtitle([num2str(1/f_p),'s'],'fontsize',18,'fontweight','bold');

save2pdf([figpath,'fk_array_response_',num2str(1/f_p),'s.pdf'],1,250);


