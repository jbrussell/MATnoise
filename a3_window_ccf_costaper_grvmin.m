% JBR - Plot the cross-spectra in the time domain for the individual
% station pairs, comparing two components. Uses same structure as
% plot_ccf_record.m
%
% NJA, 04/08/2016
% JBR, 6/17/2016 - window in time domain;
% JBR, 2/28/2018 - use less harsh cosine taper for window
% JBR, 7/8/2017 - This version only defines a minimum group velocity

clear; close all;
setup_parameters;
IsFigure = 1;

%======================= PARAMETERS =======================%comps = {'TT'}; %{'RR' 'TT' 'ZZ'}; %'ZZ'; %'RR'; %'TT'; % choose 2 components to plot

% % LOVE
% comps = {'TT'}; %{'RR' 'TT' 'ZZ'}; %'ZZ'; %'RR'; %'TT'; % choose 2 components to plot
% coperiod = [3 9]; %[2 5]; %[5 25]; %[2 8]; %[2 16]; %[5 25]; %[ 8 30 ]; % Periods to filter between
% windir = 'window3hr'; %'window3hr_LH_Zcorr'; %'window3hr_LH_whitesm'; %'window3hr_LH_Zcorr';
% % Mode Branches
% max_grv = inf; %5.5;
% min_grv = 1.4; %1.6; %2.2;

comps = {'RR'}; %{'RR' 'TT' 'ZZ'}; %'ZZ'; %'RR'; %'TT'; % choose 2 components to plot
coperiod = [3 9]; %[2 5]; %[5 25]; %[2 8]; %[2 16]; %[5 25]; %[ 8 30 ]; % Periods to filter between
windir = 'window3hr'; %'window3hr_LH_Zcorr'; %'window3hr_LH_whitesm'; %'window3hr_LH_Zcorr';
% Mode Branches
max_grv = inf; %5.5;
min_grv = 1.4; %1.6; %2.2;


IsVelLines = 1;
% WATER
H20grv = 1.4;
%==========================================================%
%xlims = [-250 250];
xlims = [-500 500];

stalist = parameters.stalist;
nsta = parameters.nsta;
winlength = parameters.winlength;
figpath = parameters.figpath;
%fig_winlength_path = [figpath,'window',num2str(winlength),'hr/fullStack/'];
% custom directory names
    fig_winlength_path = [figpath,windir,'/fullStack/'];

%------------ PATH INFORMATION -------------%
% OLD CCF
ccf_path = parameters.ccfpath;
%ccf_winlength_path = [ccf_path,'window',num2str(winlength),'hr/'];
    ccf_winlength_path = [ccf_path,windir,'/'];
ccf_singlestack_path = [ccf_winlength_path,'single/'];
ccf_daystack_path = [ccf_winlength_path,'dayStack/'];
ccf_monthstack_path = [ccf_winlength_path,'monthStack/'];
ccf_fullstack_path = [ccf_winlength_path,'fullStack/'];

ccf_stack_path = ccf_fullstack_path;

figpath = [fig_winlength_path,num2str(coperiod(1)),'_',num2str(coperiod(2)),'s/'];
% create figure directory
if ~exist(fig_winlength_path)
    mkdir(fig_winlength_path)
end
if ~exist(figpath)
    mkdir(figpath)
end

%------------ LOAD DATA AND PLOT IN TIME DOMAIN -------------%
for icomp = 1:length(comps) % loop over components
    comp = comps{icomp};
    ccf_path = [ccf_stack_path,'ccf',comp,'/',];
    npairall = 0;
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

            filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);

            if ~exist(filename,'file') % check that ccf file exists
                disp(['not exist ',filename])
                continue;
            end
            nstapair = nstapair + 1;

            %----------- LOAD DATA -------------%
            data = load(filename);
            ccf = data.coh_sum./data.coh_num;

            %----------- Frequency ==> Time domain -------------%
            N = length(ccf);
            ccf_ifft_full = ifft(ccf,N);
            
            % Distance between sta1 and sta2
            sta1sta2_dist{ista1}(nstapair) = deg2km(distance(data.stapairsinfo.lats(1),data.stapairsinfo.lons(1),data.stapairsinfo.lats(2),data.stapairsinfo.lons(2)));
            
            % Window CCF                      
            time = [-floor(N/2):1:floor(N/2)];
            ccf_ifft_full = fftshift(ccf_ifft_full); % Rearrange for windowing
            if sta1sta2_dist{ista1}(nstapair) >= 100  
                t_pos = sta1sta2_dist{ista1}(nstapair)/min_grv;
            elseif sta1sta2_dist{ista1}(nstapair) < 100
                t_pos = 100/min_grv;
            end
            Iset0_pos = (time >= -t_pos) & (time <= t_pos);
            ccf_ifft_full_pos = ccf_ifft_full;
            ccf_ifft_full_pos(Iset0_pos)= cos_taper(ccf_ifft_full_pos(Iset0_pos));
            ccf_ifft_full_pos(~Iset0_pos)= 0;

            ccf_ifft_full = ifftshift(ccf_ifft_full_pos); % Shift window back
            
            ccf_win = fft(ccf_ifft_full);
            
            %ccf_ifft = real(ifft(2*ccf([1:N/2+1]),N)); % inverse FFT to get time domain
            ccf_ifft = real(ifft(2*ccf_win([1:N/2+1]),N)); % inverse FFT to get time domain
          
            
            %rearrange and keep values corresponding to lags: -(len-1):+(len-1)
            ccf_ifft = [ccf_ifft(end-N+2:end) ; ccf_ifft(1:N)];
            
            %----------- SAVE WINDOWED CCF -------------%
            coh_sum_win = ccf_win.*data.coh_num;
            coh_sum = data.coh_sum;
            coh_num = data.coh_num;
            stapairsinfo = data.stapairsinfo;
            ccfcomp_fullstack_path = [ccf_path];                
            save(sprintf('%s%s/%s_%s_f.mat',ccfcomp_fullstack_path,sta1,sta1,sta2),'coh_sum','coh_num','stapairsinfo','coh_sum_win','max_grv','min_grv'); 

            %----------- FILTER DATA -------------%
            f1 = 1/coperiod(2);
            f2 = 1/coperiod(1);
            [b, a] = butter(2,[f1 f2]); % Butterworth Filter
            ccf_filt{ista1}{icomp}{nstapair} =  filtfilt(b,a,ccf_ifft);

            %----------- NORMALIZE CCF FUNCTION -------------%
            ccf_filt{ista1}{icomp}{nstapair} = ccf_filt{ista1}{icomp}{nstapair}/max(abs(ccf_filt{ista1}{icomp}{nstapair}));


            % Check if reverse station pair has already been plotted
            stapairinv = [sta2,'_',sta1];
            if exist('existpair','var')
                if find(strncmp(stapairinv,existpair,length(stapairinv)))
                    continue
                end
            end

            % Update some other useful variables
            dumsta2{nstapair} = sta2;
            npairall = npairall + 1; % number of total station pairs
            ccf_all{icomp}{npairall} = ccf_filt{ista1}{icomp}{nstapair} ; % cell containing all ccf
            sta1sta2_dist_all(npairall) = sta1sta2_dist{ista1}(nstapair); % vector containing distance between each station pair
            existpair(npairall) = {[sta1,'_',sta2]};


        end % ista2
    end % ista1
end % icomp

%%
clr = [0 0 0; 1 0 0; 0 0 1]; % colors for plotting [R G B]; 
% for ista1=1:nsta % loop over all stations
%     sta1=char(stalist(ista1,:));
%     if IsFigure
%         f101 = figure(101); clf; hold on;
%         for icomp = 1:length(comps) % loop over components
%             %----------- PLOT CCFs IN DISTANCE-TIME -------------%
%             for ista2 = 1: nsta % loop over station pairs
%                 N= length(ccf_ifft);
%                 time = [-N/2:N/2]; % build lagtime vector for plotting
%                 amp = 1e1; % amplitude to plot data
%                 indtime = find(abs(time)<=500); % Time index -500 to 500 seconds
%                 set(gca,'YDir','reverse');
%                 for istapair = 1: nstapair % loop over station pairs
%                     % Normalize using the surface wave amplitude
%                     tnorm = sta1sta2_dist{ista1}(istapair)/2;
%                     indtnorm = find(abs(time)<=tnorm);
%                     ccf_filt_norm = ccf_filt{ista1}{icomp}{istapair}/max(abs(ccf_filt{ista1}{icomp}{istapair}(indtnorm(1):indtnorm(end))));
%                     
%                     % Plot the normalized traces
%                     ccf_waveform = ccf_filt_norm(indtime(1):indtime(end)); % ccf at -500 to 500 seconds
%                     h1(icomp) = plot(time(indtime(1):indtime(end)),ccf_waveform*amp+sta1sta2_dist{ista1}(istapair),'-','color',clr(icomp,:)); hold on;
%                     %             text(0,stapairdist(istapair),dumsta2{istapair})
%                     %             return
%                     %             pause
%                 end
%                 xlim(xlims)
%                 xlabel('lag time (s)','fontsize',18,'fontweight','bold');
%                 ylabel('Distance (km)','fontsize',18,'fontweight','bold');
%                 title(['reference station:',sta1,'  filtered at ',num2str(coperiod(1)), ' - ',num2str(coperiod(2)),'(s)'],'fontsize',18,'fontweight','bold');
%                 set(gca,'fontsize',15);
%                 drawnow;
%             end
%         end
%         legend(h1,{comps{1}(1) comps{2}(1) comps{3}(1)},'location','northeast','fontsize',12);
%     end
% end


%%
%----------- PLOT ALL CCFs STATION PAIRS IN DISTANCE-TIME -------------%
N= length(ccf_ifft);
time = [-N/2:N/2];
amp = 1e1;
indtime = find(abs(time)<=500);

f102 = figure(102);
clf
hold on;
set(gca,'YDir','reverse');
for icomp = 1:length(comps) % loop over components
    for istapair = 1: npairall
        % Normalize using the surface wave amplitude
        ccf_filt_norm = ccf_all{icomp}{istapair}/max(abs(ccf_all{icomp}{istapair}(:)));
        
        % Plot the normalized traces
        ccf_waveform_all = ccf_filt_norm(indtime(1):indtime(end));
        h2(icomp) = plot(time(indtime(1):indtime(end)),ccf_waveform_all*amp+sta1sta2_dist_all(istapair),'-','color',clr(icomp,:)); hold on;
    end
end
% Plot lines of velocity
%plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/maxgrv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'-g');
%plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/-maxgrv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'-g');

% Plot Velocities
if IsVelLines
    % Branches
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/max_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[0 .9 0],'linewidth',2);
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/-max_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[0 .9 0],'linewidth',2);
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/min_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[255 128 0]/255,'linewidth',2);
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/-min_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[255 128 0]/255,'linewidth',2);
    
    % H20
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/H20grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[.5 .5 1],'linewidth',2);
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/-H20grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[.5 .5 1],'linewidth',2);
end
xlim(xlims)
xlabel('lag time (s)','fontsize',18,'fontweight','bold');
ylabel('Distance (km)','fontsize',18,'fontweight','bold');
title(['All non-repeated pairs, filtered at ',num2str(coperiod(1)), ' -',num2str(coperiod(2)),'(s)'],'fontsize',18,'fontweight','bold');
set(gca,'fontsize',15);
% legend(h2,{comps{1}(1) comps{2}(1) comps{3}(1)},'location','northeast','fontsize',12);

%pause;

%print(f102,'-dpdf',[figpath,'all_ccf',comps{1}(1),comps{2}(1),'2.pdf']); % Save figure
