% Plot the cross-spectra in the time domain for the individual station pairs. 
% Filter is built using a Tukey taper with sharpness controlled by costap_wid. 
% The typical butterworth filter is not precise enough for the higest frequencies
%
% https://github.com/jbrussell

clear;
if ~exist('setup_parameters_MATnoise.m')
    !cp ../setup_parameters.m ./setup_parameters_MATnoise.m
end
setup_parameters_MATnoise;
setup_parameters;
IsFigure = 0;
IsFigure_GAUS = 0; % Plot frequency domain filtered and unfiltered

%======================= PARAMETERS =======================%
comps = {'ZZ'};  % {'ZZ','RR','PP'}; 'PZ'; 'PP'; 'ZZ'; 'RR'; 'TT';
amp = 8e0;
% windir = 'window3hr';
% windir = 'window3hr_Zcorr_tiltcomp';

% Define group velocity window for calculating signal-to-noise ratio
max_grv = inf; %5.5;
min_grv = 0.7; %1.6
snr_thresh = 5; % Signal-to-noise tolerance

% Define time axis limits for cutting waveforms
xlims = [0 500];
ylims = [0 450];


%%% --- Parameters to build up gaussian filters --- %%% 
% (effects the width of the filter in the frequency domain)
costap_wid = 0.2; % 0 => box filter; 1 => Hann window

isplotwin = 1; %1;
isploth20 = 0;
isfigure_snr = 0;

h20_grv = 1.5;
%==========================================================%

dt = parameters.dt;
stalist = parameters.stalist;
nsta = parameters.nsta;
nsta = length(stalist);
winlength = parameters.winlength;
figpath = parameters.figpath;
%fig_winlength_path = [figpath,'window',num2str(winlength),'hr/fullStack/'];
% custom directory names
    fig_winlength_path = [figpath,windir,'/fullStack/'];

%------------ PATH INFORMATION -------------%
% ccf_path = parameters.ccfpath;
%ccf_winlength_path = [ccf_path,'window',num2str(winlength),'hr/'];
    ccf_winlength_path = [ccfpath,windir,'/'];
ccf_singlestack_path = [ccf_winlength_path,'single/'];
ccf_daystack_path = [ccf_winlength_path,'dayStack/'];
ccf_monthstack_path = [ccf_winlength_path,'monthStack/'];
ccf_fullstack_path = [ccf_winlength_path,'fullStack/'];

ccf_stack_path = ccf_fullstack_path;

% create figure directory
if ~exist(fig_winlength_path)
    mkdir(fig_winlength_path)
end
if ~exist(figpath)
    mkdir(figpath)
end

%% Load Depths
STAS = stalist;
LATS = stalat;
LONS = stalon;
DEPTHS = staz;

M = [];
Delta = [];
t = [];

for icomp = 1:length(comps)
    comp = comps{icomp};
    %%
    ccf_path = [ccf_stack_path,'ccf',comp,'/',];
    npairall = 0;
    %------------ LOAD DATA AND PLOT IN TIME DOMAIN -------------%
    for ista1=1:nsta % loop over all stations
        sta1=char(stalist(ista1,:));
        sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
        nstapair = 0;
        for ista2 = 1: nsta % loop over station pairs
            sta2 = char(stalist(ista2,:));

            % if same station, skip
            if(strcmp(sta1,sta2))
                continue
            end

            % Check if reverse station pair has already been plotted
            npairall = npairall+1;
            stapairinv = [sta2,'_',sta1];
            if exist('existpair','var')
                if find(strncmp(stapairinv,existpair,length(stapairinv)))
                    continue
                end
            end
            existpair(npairall) = {[sta1,'_',sta2]};

            filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);

            if ~exist(filename,'file') % check that ccf file exists
                disp(['not exist ',filename])
                continue;
            end
            nstapair = nstapair + 1;

            %----------- LOAD DATA -------------%
            data = load(filename);
            ccf = data.coh_sum./data.coh_num;
            ccf(isnan(ccf)) = 0;
            if size(ccf,1) == 1
                ccf = ccf';
            end
            
            % Calculate SNR
            r = distance(data.stapairsinfo.lats(1),data.stapairsinfo.lons(1),data.stapairsinfo.lats(2),data.stapairsinfo.lons(2),referenceEllipsoid('GRS80'))/1000;
            [snr, signal_ind] = calc_SNR(ccf,min_grv,max_grv,r,isfigure_snr);
            if snr < snr_thresh
                continue
            end

            %%
            %----------- Frequency ==> Time domain -------------%
            N = length(ccf);
            ccf_ifft = real(ifft(ccf,N)); % inverse FFT to get time domain
            ccf_ifft = fftshift(ccf_ifft); % rearrange values as [-lag lag]
            ccf_ifft = detrend(ccf_ifft);
            ccf_ifft = cos_taper(ccf_ifft);

            % Build time axis [-lag:0:lag] and index positive and negative waveforms
            time = ([0:N-1]-floor(N/2))*dt;  % build lagtime vector for plotting
            time = [time(time<0), time(time>=0)];
            indtime_pos = find(abs(time)<=xlims(2) & time>=0);
            indtime_neg = find(abs(time)<=xlims(2) & time<=0);

            ccf_ifft_pos = (detrend(ccf_ifft(indtime_pos)));
            ccf_ifft_neg = flip((detrend(ccf_ifft(indtime_neg))));            
%             ccf_ifft_pos = cos_taper(detrend(ccf_ifft(indtime_pos)));
%             ccf_ifft_neg = flip(cos_taper(detrend(ccf_ifft(indtime_neg))));

    %         % flip -lag to positive +lag
    %         M = [M; ccf_ifft_pos; ccf_ifft_neg ];
    %         Delta = [Delta; r; r ];
    %         t = time(indtime_pos);

            % average +lag and -lag
            M = [M; mean([ccf_ifft_pos; ccf_ifft_neg],1)];
            Delta = [Delta; r ];
            t = time(indtime_pos);

    %         % Keep -lag:+lag
    %         M = [M; cos_taper(detrend(ccf_ifft(abs(time)<=xlims(2)))) ];
    %         Delta = [Delta; r ];
    %         t = time(abs(time)<=xlims(2));


        end % ista2
    end % ista1
end

% Sort data by distance
[Delta, I_srt] = sort(Delta);
M = M(I_srt,:);

% Apply group velocity window
I_keep = t<=Delta./min_grv & t>=Delta./max_grv;
win_mat = zeros(size(I_keep));
for itr = 1:length(Delta)
    I_onesrow = I_keep(itr,:)==1;
    win_row = I_keep(itr,I_onesrow);
    win_mat(itr,I_onesrow) = cos_taper(win_row);
end
M_win = M.*win_mat;

%% %----------- PLOT ALL CCFs STATION PAIRS IN DISTANCE-TIME -------------%
f102 = figure(102);
set(gcf, 'Color', 'w');
clf
hold on;
set(gca,'YDir','reverse');
plot(t,M./max(abs(M),[],2)*amp+Delta,'-k','linewidth',1);
plot(t,M_win./max(abs(M),[],2)*amp+Delta,'-r','linewidth',0.5);
xlim([min(t) max(t)])
% xlim([0 max(xlims)])
ylim(ylims);
xlabel('lag time (s)','fontsize',18);
ylabel('Distance (km)','fontsize',18);
set(gca,'fontsize',15);

% Plot Velocities
if isplotwin
    % Branches
    plot([min(Delta) max(Delta)]/max_grv,[min(Delta) max(Delta)],'color',[1 0 0],'linewidth',2);
    plot([min(Delta) max(Delta)]/-max_grv,[min(Delta) max(Delta)],'color',[1 0 0],'linewidth',2);
    plot([min(Delta) max(Delta)]/min_grv,[min(Delta) max(Delta)],'color',[1 0 0],'linewidth',2);
    plot([min(Delta) max(Delta)]/-min_grv,[min(Delta) max(Delta)],'color',[1 0 0],'linewidth',2);    
end


odata = [datapath,'noise_',comps{:},'.mat'];
if ~exist(datapath)
    mkdir(datapath)
end
save(odata,'M','M_win','max_grv','min_grv','Delta','t');