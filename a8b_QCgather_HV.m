% Gather H/V ellipticity values for each station and do some statistics to
% determine station averages. Must first run a8a_calc_HV_amp_ftan.
%
% jbrussell - 10/2025
clear all; close all;
setup_parameters;
addpath('./ray_tomo/tomo_functions/')

IsFigure = 1;
isoverwrite = 1; % overwrite results?
isoutput = 1; % save output?

%% Parameters from FTAN script
windir = 'window3hr_ellip';
frange_fit = [1/25 1/4]; % Frequency range to estimate grv
opts.nBranches = 1; %    (default 2)     % # dispersion branches to pick

%QC Parameters
snr_thresh = 8; % minimum SNR to consider. If any component (ZZ,RR,ZR,RZ) is less than this value then measurement is dropped
wl_thresh = 1; % minimum wavelength to consider
dphi_thresh = 90+[-45 +45]; % ensure Rayleigh wave behavior
tg_period_frac_thresh = 5; %1/4; % ensure traveltime of RR,RZ,ZR within fraction of period relative to ZZ
hv_RZ_diff_thresh = 0.1; % (default 10%) fractional difference between H/V estimates for vertical and radial forces
is_3sig_removal = 1; % (1:yes | 0:no) Removal outliers greater than 3 standard deviations from mean

%%
% dt = parameters.dt;
stalist = parameters.stalist;
nsta = parameters.nsta;
winlength = parameters.winlength;
figpath = parameters.figpath;
%fig_winlength_path = [figpath,'window',num2str(winlength),'hr/fullStack/'];
% custom directory names
fig_winlength_path = [figpath,windir,'/fullStack/'];

%% Setup directories

HV_in_path = ['./HV_ellip/',windir,'/fullStack/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
HV_out_path = ['./HV_ellip_stats/',windir,'/fullStack/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
if ~exist(HV_out_path)
    mkdir(HV_out_path)
end

% figure output path
hv_fig_path = ['./figs/',windir,'/fullStack/HV_ellip_stats/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
if ~exist(hv_fig_path)
    mkdir(hv_fig_path);
end

%% Loop over stations and calculate ellipticity

for ista=1:nsta % loop over all stations
    % Get station of interest
    sta=char(stalist(ista,:));

    HV_out_filename = sprintf('%s/%s_hvstats.mat',HV_out_path,sta);
    if ~isoverwrite
        if exist(HV_out_filename)
            disp(['Already processed ',sta,' ... skipping']);
            continue;
        end
    end

    % Initialize arrays
    hv_vals = [];
    hv_RZdiff_vals = [];
    tg_difffromZZ_vals = [];
    hv_stas = {};
    r_vals = [];
    wl_vals = [];
    snr_vals = [];
    dphi_vals = [];

    % Start with files in which station of interest is sta1
    files = dir([HV_in_path,'/',sta,'_*_hv.mat']);
    for ii = 1:length(files)
        temp = load([HV_in_path,'/',files(ii).name]);
        hv = temp.hv;
        periods = hv.periods; 
        sta1 = hv.stapairsinfo.stanames{1};
        sta2 = hv.stapairsinfo.stanames{2};
        hv_1_stack = hv.sta1.RZ_ZZ_stack;
        hv_2_stack = hv.sta1.RR_ZR_stack;
        hv_av_stack = 0.5 * (hv_1_stack + hv_2_stack);
        hv_RZdiff_stack = 1 - hv_1_stack./hv_2_stack;
        tg_RR_difffromZZ_stack = abs(hv.stapairsinfo.r./hv.grv.RR_stack - hv.stapairsinfo.r./hv.grv.ZZ_stack);
        tg_RZ_difffromZZ_stack = abs(hv.stapairsinfo.r./hv.grv.RZ_stack - hv.stapairsinfo.r./hv.grv.ZZ_stack);
        tg_ZR_difffromZZ_stack = abs(hv.stapairsinfo.r./hv.grv.ZR_stack - hv.stapairsinfo.r./hv.grv.ZZ_stack);
        hv_vals = [hv_vals; hv_av_stack(:)'];
        hv_RZdiff_vals = [hv_RZdiff_vals; hv_RZdiff_stack(:)'];
        tg_difffromZZ_vals = [tg_difffromZZ_vals; max([tg_RR_difffromZZ_stack(:)'; tg_RZ_difffromZZ_stack(:)'; tg_ZR_difffromZZ_stack(:)'])];
        hv_stas = [hv_stas; [sta1,'-',sta2]];
        r_vals = [r_vals; hv.stapairsinfo.r];
        wl_vals = [wl_vals; hv.stapairsinfo.r./(hv.grv.ZZ_stack(:)'.*periods(:)')];
%         snr_vals = [snr_vals; 0.25*(hv.snr.ZZ_stack + hv.snr.RR_stack + hv.snr.ZR_stack + hv.snr.RZ_stack)]; % average SNR
        snr_vals = [snr_vals; min([hv.snr.ZZ_stack; hv.snr.RR_stack; hv.snr.ZR_stack; hv.snr.RZ_stack])]; % minimum SNR
        dphi_vals = [dphi_vals; angmean(pi/180*[hv.sta1.dphi_RZ_ZZ_stack(:)'; hv.sta1.dphi_RR_ZR_stack(:)'])*180/pi];
        
    end 

    % Now loop over files in which station of interest is sta2
    files = dir([HV_in_path,'/','*_',sta,'_hv.mat']);
    for ii = 1:length(files)
        temp = load([HV_in_path,'/',files(ii).name]);
        hv = temp.hv;
        sta1 = hv.stapairsinfo.stanames{1};
        sta2 = hv.stapairsinfo.stanames{2};
        hv_1_stack = hv.sta2.ZR_ZZ_stack;
        hv_2_stack = hv.sta2.RR_RZ_stack;
        hv_av_stack = 0.5 * (hv_1_stack + hv_2_stack);
        hv_RZdiff_stack = 1 - hv_1_stack./hv_2_stack;
        tg_RR_difffromZZ_stack = abs(hv.stapairsinfo.r./hv.grv.RR_stack - hv.stapairsinfo.r./hv.grv.ZZ_stack);
        tg_RZ_difffromZZ_stack = abs(hv.stapairsinfo.r./hv.grv.RZ_stack - hv.stapairsinfo.r./hv.grv.ZZ_stack);
        tg_ZR_difffromZZ_stack = abs(hv.stapairsinfo.r./hv.grv.ZR_stack - hv.stapairsinfo.r./hv.grv.ZZ_stack);
        hv_vals = [hv_vals; hv_av_stack(:)'];
        hv_RZdiff_vals = [hv_RZdiff_vals; hv_RZdiff_stack(:)'];
        tg_difffromZZ_vals = [tg_difffromZZ_vals; max([tg_RR_difffromZZ_stack(:)'; tg_RZ_difffromZZ_stack(:)'; tg_ZR_difffromZZ_stack(:)'])];
        hv_stas = [hv_stas; [sta1,'-',sta2]];
        r_vals = [r_vals; hv.stapairsinfo.r];
        wl_vals = [wl_vals; hv.stapairsinfo.r./(hv.grv.ZZ_stack(:)'.*periods(:)')];
%         snr_vals = [snr_vals; 0.25*(hv.snr.ZZ_stack + hv.snr.RR_stack + hv.snr.ZR_stack + hv.snr.RZ_stack)]; % average SNR
        snr_vals = [snr_vals; min([hv.snr.ZZ_stack; hv.snr.RR_stack; hv.snr.ZR_stack; hv.snr.RZ_stack])]; % minimum SNR
        dphi_vals = [dphi_vals; angmean(pi/180*[hv.sta2.dphi_ZR_ZZ_stack(:)'; hv.sta2.dphi_RR_RZ_stack(:)'])*180/pi];
    end

    if isempty(hv_vals)
        disp(['No measurements for ',sta,'... skipping'])
        continue
    end

    %% Do quality control
    I_QCpass = snr_vals >= snr_thresh & ...
               wl_vals >= wl_thresh & ...
               (abs(dphi_vals)>=dphi_thresh(1) & abs(dphi_vals)<=dphi_thresh(2)) & ...
               abs(tg_difffromZZ_vals)./repmat(periods,size(tg_difffromZZ_vals,1),1) <= tg_period_frac_thresh & ...
               hv_RZdiff_vals <= hv_RZ_diff_thresh;
    
    hv_vals(~I_QCpass) = nan;
    dphi_vals(~I_QCpass) = nan;

    % Remove large outliers (> 3std)
    if is_3sig_removal
        Ioutliers = abs(hv_vals-repmat(nanmean(hv_vals,1),size(hv_vals,1),1)) >= 3*nanstd(hv_vals,[],1);
    end
    
    Igood = I_QCpass & ~Ioutliers;
    hv_vals(~Igood) = nan;
    dphi_vals(~Igood) = nan;

    %% Get and save stats
    hv_stats.hv_med = nanmedian(hv_vals,1);
    hv_stats.hv_mean = nanmean(hv_vals,1);
    hv_stats.hv_std = nanstd(hv_vals,[],1);
    for iper = 1:length(periods)
        hv_stats.dphi_mean(1,iper) = angmean(pi/180*dphi_vals(Igood(:,iper),iper))*180/pi;
        hv_stats.dphi_std(1,iper) = angstd(pi/180*dphi_vals(Igood(:,iper),iper))*180/pi;
    end
    hv_stats.periods = periods;
    hv_stats.hv_vals = hv_vals;
    hv_stats.hv_stas = hv_stas;
    hv_stats.Igood = Igood;
    hv_stats.r_vals = r_vals;
    hv_stats.wl_vals = wl_vals;
    hv_stats.snr_vals = snr_vals;
    hv_stats.dphi_vals = dphi_vals;
    hv_stats.QC_params.snr_thresh = snr_thresh;
    hv_stats.QC_params.wl_thresh = wl_thresh;
    hv_stats.QC_params.dphi_thresh = dphi_thresh;

    if isoutput
        save(HV_out_filename,'hv_stats');
    end
        
    %% Plot station results
    if IsFigure
        figure(10000); clf;
        set(gcf,'position',[137  226  978  724],'color','w')
        
        subplot(3,1,1);
        box on; hold on;
        errorbar(hv_stats.periods,hv_stats.hv_med,hv_stats.hv_std,'-or','linewidth',2)
        xlabel('Period (s)');
        ylabel('H/V');
        title([sta])
%         legend('location','eastoutside')
        set(gca,'fontsize',15,'linewidth',1.5);
        xlim([min(hv_stats.periods) max(hv_stats.periods)]);

        subplot(3,1,2);
        box on; hold on;
        errorbar(hv_stats.periods,hv_stats.dphi_mean,hv_stats.dphi_std,'-or','linewidth',2)
        yline(90,'--k');
        yline(-90,'--k');
        xlabel('Period (s)');
        ylabel('Phase difference (Z-R)');
        title([sta])
%         legend('location','eastoutside')
        set(gca,'fontsize',15,'linewidth',1.5);
        xlim([min(hv_stats.periods) max(hv_stats.periods)]);
        ylim([-180 180]);
        yticks([-180:90:180]);
        
        Nplot = 5;
        ipers = round(length(hv_stats.periods)/Nplot * [1:Nplot]);
        for ii = 1:length(ipers)
            subplot(3,Nplot,2*Nplot+ii);
            box on; hold on;
            histogram(hv_stats.hv_vals(:,ipers(ii)),10);
            xlabel('H/V')
            title([num2str(hv_stats.periods(ipers(ii))),'s'])
            set(gca,'fontsize',15,'linewidth',1.5);
        end
        if isoutput
            save2pdf([hv_fig_path,'/',sta,'_HVstats.pdf'],10000,300);
        end
    end
        
end