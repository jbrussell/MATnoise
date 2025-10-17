% Use the outputs of FTAN for ZZ, RR, ZR, and RZ for estimation of
% Rayleigh-wave ellipticity following Lin et al. (2014). This script reads
% in the FTAN mat structure from a6_b_measure_ftan_grv and requires it to
% have been run for ZZ, RR, ZR, and RZ. Ellipticity is measured for the
% positive (causal), negative (acausal), and stacked ccf results.
%
% The instantaneous phase difference between Z-R components is also calculated
% from the narrow-band filtered waveforms. Instantaneous phase is extracted 
% at the group travel time estimated from FTAN. The convention used here 
% is as follows:
%           -90º = retrograde (R leads Z)
%           +90º = prograde (Z leads R)
% Because the ambient noise cross-correlation assumes R points from sta1 to
% sta2, when we estimate phase difference for ellipticity measured at sta1,
% the horizontal component must be multiplied by -1.
%
% jbrussell - 10/2025
clear all; close all;
setup_parameters;

IsFigure = 1;
isoverwrite = 1; % overwrite results?
isoutput = 1; % save output?

%% Parameters from FTAN script
windir = 'window3hr_ellip';
frange_fit = [1/25 1/4]; % Frequency range to estimate grv
opts.nBranches = 1; %    (default 2)     % # dispersion branches to pick

% FTAN directory names (should not need to edit these)
path2ftan_ZZ = ['./ftan/',windir,'/fullStack/ftanZZ/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
path2ftan_RR = ['./ftan/',windir,'/fullStack/ftanRR/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
path2ftan_RZ = ['./ftan/',windir,'/fullStack/ftanRZ/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
path2ftan_ZR = ['./ftan/',windir,'/fullStack/ftanZR/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];


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

HV_out_path = ['./HV_ellip/',windir,'/fullStack/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
if ~exist(HV_out_path)
    mkdir(HV_out_path)
end

% figure output path
hv_fig_path = ['./figs/',windir,'/fullStack/HV_ellip/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
if ~exist(hv_fig_path)
    mkdir(hv_fig_path);
end

%% Loop over stations and calculate ellipticity

for ista1=1:nsta % loop over all stations
    sta1=char(stalist(ista1,:));
    
    for ista2 = 1: nsta % loop over station pairs
        sta2 = char(stalist(ista2,:));
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        
        HV_out_filename = sprintf('%s/%s_%s_hv.mat',HV_out_path,sta1,sta2);
        if ~isoverwrite
            if exist(HV_out_filename)
                disp(['Already processed ',sta1,'-',sta2,' ... skipping']);
                continue;
            end
        end
                
        filename_ZZ = sprintf('%s/%s_%s_ftan.mat',path2ftan_ZZ,sta1,sta2);
        filename_RR = sprintf('%s/%s_%s_ftan.mat',path2ftan_RR,sta1,sta2);
        filename_RZ = sprintf('%s/%s_%s_ftan.mat',path2ftan_RZ,sta1,sta2);
        filename_ZR = sprintf('%s/%s_%s_ftan.mat',path2ftan_ZR,sta1,sta2);
        
        if ~exist(filename_ZZ,'file') % check that ftan files exist
            disp(['not exist ',filename_ZZ])
            continue;
        end
        if ~exist(filename_RR,'file') % check that ftan files exist
            disp(['not exist ',filename_RR])
            continue;
        end
        if ~exist(filename_RZ,'file') % check that ftan files exist
            disp(['not exist ',filename_RZ])
            continue;
        end
        if ~exist(filename_ZR,'file') % check that ftan files exist
            disp(['not exist ',filename_ZR])
            continue;
        end

        temp = load(filename_ZZ);
        ftan_ZZ = temp.ftan;
        temp = load(filename_RR);
        ftan_RR = temp.ftan;
        temp = load(filename_RZ);
        ftan_RZ = temp.ftan;
        temp = load(filename_ZR);
        ftan_ZR = temp.ftan; clear temp;

        periods = ftan_ZZ.periods;

        %% Do positive (causal) side first
        amp_ZZ_pos = ftan_ZZ.pos.amp;
        amp_RR_pos = ftan_RR.pos.amp;
        amp_RZ_pos = ftan_RZ.pos.amp;
        amp_ZR_pos = ftan_ZR.pos.amp;
        time = ftan_RZ.time(ftan_RZ.time>=0);
        
        hv.sta1.RZ_ZZ_pos = amp_RZ_pos ./ amp_ZZ_pos;
        hv.sta1.RR_ZR_pos = amp_RR_pos ./ amp_ZR_pos;
        hv.sta2.ZR_ZZ_pos = amp_ZR_pos ./ amp_ZZ_pos;
        hv.sta2.RR_RZ_pos = amp_RR_pos ./ amp_RZ_pos;

        % Get phase relationships
        % RZ/ZZ (multiply R by -1 because of coordinate convention: R points sta1-->sta2 by default)
        dphi = calc_inst_phase_diff(ftan_ZZ.pos.ccf,-1*ftan_RZ.pos.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_ZZ.pos.tg + ftan_RZ.pos.tg); % average group travel time
        hv.sta1.dphi_RZ_ZZ_pos = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta1.dphi_RZ_ZZ_pos(iper) = dphi(itg,iper);
        end
        % RR/ZR (multiply R by -1 because of coordinate convention: R points sta1-->sta2 by default)
        dphi = calc_inst_phase_diff(ftan_ZR.pos.ccf,-1*ftan_RR.pos.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_ZR.pos.tg + ftan_RR.pos.tg); % average group travel time
        hv.sta1.dphi_RR_ZR_pos = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta1.dphi_RR_ZR_pos(iper) = dphi(itg,iper);
        end
        % ZR/ZZ
        dphi = calc_inst_phase_diff(ftan_ZZ.pos.ccf,ftan_ZR.pos.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_ZZ.pos.tg + ftan_ZR.pos.tg); % average group travel time
        hv.sta2.dphi_ZR_ZZ_pos = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta2.dphi_ZR_ZZ_pos(iper) = dphi(itg,iper);
        end
        % RR/RZ
        dphi = calc_inst_phase_diff(ftan_RZ.pos.ccf,ftan_RR.pos.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_RZ.pos.tg + ftan_RR.pos.tg); % average group travel time
        hv.sta2.dphi_RR_RZ_pos = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta2.dphi_RR_RZ_pos(iper) = dphi(itg,iper);
        end
        
        %% Now do negative (acausal) side
        amp_ZZ_neg = ftan_ZZ.neg.amp;
        amp_RR_neg = ftan_RR.neg.amp;
        amp_RZ_neg = ftan_RZ.neg.amp;
        amp_ZR_neg = ftan_ZR.neg.amp;
        
        hv.sta1.RZ_ZZ_neg = amp_RZ_neg ./ amp_ZZ_neg;
        hv.sta1.RR_ZR_neg = amp_RR_neg ./ amp_ZR_neg;
        hv.sta2.ZR_ZZ_neg = amp_ZR_neg ./ amp_ZZ_neg;
        hv.sta2.RR_RZ_neg = amp_RR_neg ./ amp_RZ_neg;

        % Get phase relationships
        % RZ/ZZ (multiply R by -1 because of coordinate convention: R points sta1-->sta2 by default)
        dphi = calc_inst_phase_diff(ftan_ZZ.neg.ccf,-1*ftan_RZ.neg.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_ZZ.neg.tg + ftan_RZ.neg.tg); % average group travel time
        hv.sta1.dphi_RZ_ZZ_neg = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta1.dphi_RZ_ZZ_neg(iper) = dphi(itg,iper);
        end
        % RR/ZR (multiply R by -1 because of coordinate convention: R points sta1-->sta2 by default)
        dphi = calc_inst_phase_diff(ftan_ZR.neg.ccf,-1*ftan_RR.neg.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_ZR.neg.tg + ftan_RR.neg.tg); % average group travel time
        hv.sta1.dphi_RR_ZR_neg = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta1.dphi_RR_ZR_neg(iper) = dphi(itg,iper);
        end
        % ZR/ZZ
        dphi = calc_inst_phase_diff(ftan_ZZ.neg.ccf,ftan_ZR.neg.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_ZZ.neg.tg + ftan_ZR.neg.tg); % average group travel time
        hv.sta2.dphi_ZR_ZZ_neg = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta2.dphi_ZR_ZZ_neg(iper) = dphi(itg,iper);
        end
        % RR/RZ
        dphi = calc_inst_phase_diff(ftan_RZ.neg.ccf,ftan_RR.neg.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_RZ.neg.tg + ftan_RR.neg.tg); % average group travel time
        hv.sta2.dphi_RR_RZ_neg = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta2.dphi_RR_RZ_neg(iper) = dphi(itg,iper);
        end

        %% Now do the average stack
        amp_ZZ_stack = ftan_ZZ.stack.amp;
        amp_RR_stack = ftan_RR.stack.amp;
        amp_RZ_stack = ftan_RZ.stack.amp;
        amp_ZR_stack = ftan_ZR.stack.amp;
        
        hv.sta1.RZ_ZZ_stack = amp_RZ_stack ./ amp_ZZ_stack;
        hv.sta1.RR_ZR_stack = amp_RR_stack ./ amp_ZR_stack;
        hv.sta2.ZR_ZZ_stack = amp_ZR_stack ./ amp_ZZ_stack;
        hv.sta2.RR_RZ_stack = amp_RR_stack ./ amp_RZ_stack;

        % Get phase relationships
        % RZ/ZZ (multiply R by -1 because of coordinate convention: R points sta1-->sta2 by default)
        dphi = calc_inst_phase_diff(ftan_ZZ.stack.ccf,-1*ftan_RZ.stack.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_ZZ.stack.tg + ftan_RZ.stack.tg); % average group travel time
        hv.sta1.dphi_RZ_ZZ_stack = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta1.dphi_RZ_ZZ_stack(iper) = dphi(itg,iper);
        end
        % RR/ZR (multiply R by -1 because of coordinate convention: R points sta1-->sta2 by default)
        dphi = calc_inst_phase_diff(ftan_ZR.stack.ccf,-1*ftan_RR.stack.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_ZR.stack.tg + ftan_RR.stack.tg); % average group travel time
        hv.sta1.dphi_RR_ZR_stack = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta1.dphi_RR_ZR_stack(iper) = dphi(itg,iper);
        end
        % ZR/ZZ
        dphi = calc_inst_phase_diff(ftan_ZZ.stack.ccf,ftan_ZR.stack.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_ZZ.stack.tg + ftan_ZR.stack.tg); % average group travel time
        hv.sta2.dphi_ZR_ZZ_stack = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta2.dphi_ZR_ZZ_stack(iper) = dphi(itg,iper);
        end
        % RR/RZ
        dphi = calc_inst_phase_diff(ftan_RZ.stack.ccf,ftan_RR.stack.ccf); % -ive = retrograde; +ive = prograde
        tg_av = 0.5*(ftan_RZ.stack.tg + ftan_RR.stack.tg); % average group travel time
        hv.sta2.dphi_RR_RZ_stack = [];
        for iper = 1:length(periods)
            [~,itg] = min(abs(time-tg_av(iper)));
            hv.sta2.dphi_RR_RZ_stack(iper) = dphi(itg,iper);
        end

        %% Save some SNR metrics
        hv.snr.ZZ_pos = ftan_ZZ.pos.snr;
        hv.snr.RR_pos = ftan_RR.pos.snr;
        hv.snr.RZ_pos = ftan_RZ.pos.snr;
        hv.snr.ZR_pos = ftan_ZR.pos.snr;
        hv.snr.ZZ_neg = ftan_ZZ.neg.snr;
        hv.snr.RR_neg = ftan_RR.neg.snr;
        hv.snr.RZ_neg = ftan_RZ.neg.snr;
        hv.snr.ZR_neg = ftan_ZR.neg.snr;
        hv.snr.ZZ_stack = ftan_ZZ.snr;
        hv.snr.RR_stack = ftan_RR.snr;
        hv.snr.RZ_stack = ftan_RZ.snr;
        hv.snr.ZR_stack = ftan_ZR.snr;

        % Save group velocities
        hv.grv.ZZ_pos = ftan_ZZ.pos.grv;
        hv.grv.RR_pos = ftan_RR.pos.grv;
        hv.grv.RZ_pos = ftan_RZ.pos.grv;
        hv.grv.ZR_pos = ftan_ZR.pos.grv;
        hv.grv.ZZ_neg = ftan_ZZ.neg.grv;
        hv.grv.RR_neg = ftan_RR.neg.grv;
        hv.grv.RZ_neg = ftan_RZ.neg.grv;
        hv.grv.ZR_neg = ftan_ZR.neg.grv;
        hv.grv.ZZ_stack = ftan_ZZ.stack.grv;
        hv.grv.RR_stack = ftan_RR.stack.grv;
        hv.grv.RZ_stack = ftan_RZ.stack.grv;
        hv.grv.ZR_stack = ftan_ZR.stack.grv;

        % Save station info
        hv.stapairsinfo = ftan_ZZ.stapairsinfo;
        hv.periods = periods;

        if isoutput
            save(HV_out_filename,'hv');
        end
        
        %% Plot station pair HV
        if IsFigure
            figure(9998); clf;
            set(gcf,'position',[115   203   938   753],'color','w')

            subplot(2,1,1);
            box on; hold on;
            h9999(1) = plot(periods,hv.sta1.RZ_ZZ_pos,'-ob','DisplayName',[sta1,' RZ/ZZ +']);
            plot(periods,hv.sta1.RZ_ZZ_neg,':ob','DisplayName',[sta1,' RZ/ZZ -']);
            plot(periods,hv.sta1.RZ_ZZ_stack,'--ob','DisplayName',[sta1,' RZ/ZZ stack']);
            plot(periods,hv.sta1.RR_ZR_pos,'-oc','DisplayName',[sta1,' RR/ZR +']);
            plot(periods,hv.sta1.RR_ZR_neg,':oc','DisplayName',[sta1,' RR/ZR -']);
            plot(periods,hv.sta1.RR_ZR_stack,'--oc','DisplayName',[sta1,' RR/ZR stack']);
            h9999(2) = plot(periods,hv.sta2.ZR_ZZ_pos,'-or','DisplayName',[sta2,' ZR/ZZ +']);
            plot(periods,hv.sta2.ZR_ZZ_neg,':or','DisplayName',[sta2,' ZR/ZZ -']);
            plot(periods,hv.sta2.ZR_ZZ_stack,'--or','DisplayName',[sta2,' ZR/ZZ stack']);
            plot(periods,hv.sta2.RR_RZ_pos,'-om','DisplayName',[sta2,' RR/RZ +']);
            plot(periods,hv.sta2.RR_RZ_neg,':om','DisplayName',[sta2,' RR/RZ -']);
            plot(periods,hv.sta2.RR_RZ_stack,'--om','DisplayName',[sta2,' RR/RZ stack']);
            xlabel('Period (s)');
            ylabel('H/V');
            xlim([min(periods) max(periods)])
            title([sta1,'-',sta2])
            legend('location','eastoutside')
            set(gca,'fontsize',15,'linewidth',1.5);

            subplot(2,1,2);
            box on; hold on;
            h9999(1) = plot(periods,hv.sta1.dphi_RZ_ZZ_pos,'-ob','DisplayName',[sta1,' RZ/ZZ +']);
            plot(periods,hv.sta1.dphi_RZ_ZZ_neg,':ob','DisplayName',[sta1,' RZ/ZZ -']);
            plot(periods,hv.sta1.dphi_RZ_ZZ_stack,'--ob','DisplayName',[sta1,' RZ/ZZ stack']);
            plot(periods,hv.sta1.dphi_RR_ZR_pos,'-oc','DisplayName',[sta1,' RR/ZR +']);
            plot(periods,hv.sta1.dphi_RR_ZR_neg,':oc','DisplayName',[sta1,' RR/ZR -']);
            plot(periods,hv.sta1.dphi_RR_ZR_stack,'--oc','DisplayName',[sta1,' RR/ZR stack']);
            h9999(2) = plot(periods,hv.sta2.dphi_ZR_ZZ_pos,'-or','DisplayName',[sta2,' ZR/ZZ +']);
            plot(periods,hv.sta2.dphi_ZR_ZZ_neg,':or','DisplayName',[sta2,' ZR/ZZ -']);
            plot(periods,hv.sta2.dphi_ZR_ZZ_stack,'--or','DisplayName',[sta2,' ZR/ZZ stack']);
            plot(periods,hv.sta2.dphi_RR_RZ_pos,'-om','DisplayName',[sta2,' RR/RZ +']);
            plot(periods,hv.sta2.dphi_RR_RZ_neg,':om','DisplayName',[sta2,' RR/RZ -']);
            plot(periods,hv.sta2.dphi_RR_RZ_stack,'--om','DisplayName',[sta2,' RR/RZ stack']);
            yline(90,'--k');
            yline(-90,'--k');
            xlabel('Period (s)');
            ylabel('Phase difference (Z-R)');
            ylim([-180 180]);
            xlim([min(periods) max(periods)])
            title([sta1,'-',sta2])
            yticks([-180:30:180]);
            legend('location','eastoutside')
            set(gca,'fontsize',15,'linewidth',1.5);

            if isoutput
                save2pdf([hv_fig_path,'/',sta1,'_',sta2,'_HVcurves.pdf'],9998,300);
            end


            %% Plot time domain

            figure(9999); clf;
            set(gcf,'position',[23   492   908   419],'color','w');
            iper = round(length(periods)/2);

            subplot(2,4,[1 2]); box on; hold on;
            plot(time,ftan_ZZ.stack.ccf(:,iper),'-r','linewidth',1.5);
            plot(time,-1*ftan_RZ.stack.ccf(:,iper),'-b','linewidth',1.5);
            legend({'ZZ';'RZ'})
            xlim([0 500])
            xlabel('Lag Time')
            title([num2str(periods(iper)),' s'])
            set(gca,'linewidth',1.5,'fontsize',12)

            subplot(2,4,[5 6]); box on; hold on;
            plot(time,ftan_ZZ.stack.ccf(:,iper),'-r','linewidth',1.5);
            plot(time,ftan_ZR.stack.ccf(:,iper),'-b','linewidth',1.5);
            legend({'ZZ';'ZR'})
            xlim([0 500])
            xlabel('Lag Time')
            title([num2str(periods(iper)),' s'])
            set(gca,'linewidth',1.5,'fontsize',12)

            subplot(2,4,3); box on;
            plot(-1*ftan_RZ.stack.ccf(:,iper),ftan_ZZ.stack.ccf(:,iper),'-r','linewidth',1.5);
            xlabel('RZ');
            ylabel('ZZ');
            axis square;
            axis equal;
            title([sta1])
            set(gca,'linewidth',1.5,'fontsize',12)

            subplot(2,4,4); box on;
            plot(-1*ftan_RR.stack.ccf(:,iper),ftan_ZR.stack.ccf(:,iper),'-r','linewidth',1.5);
            xlabel('RR');
            ylabel('ZR');
            axis square;
            axis equal;
            title([sta1])
            set(gca,'linewidth',1.5,'fontsize',12)

            subplot(2,4,7); box on;
            plot(ftan_ZR.stack.ccf(:,iper),ftan_ZZ.stack.ccf(:,iper),'-r','linewidth',1.5);
            xlabel('ZR');
            ylabel('ZZ');
            axis square;
            axis equal;
            title([sta2])
            set(gca,'linewidth',1.5,'fontsize',12)

            subplot(2,4,8); box on;
            plot(ftan_RR.stack.ccf(:,iper),ftan_RZ.stack.ccf(:,iper),'-r','linewidth',1.5);
            xlabel('RR');
            ylabel('RZ');
            axis square;
            axis equal;
            title([sta2])
            set(gca,'linewidth',1.5,'fontsize',12)

            if isoutput
                save2pdf([hv_fig_path,'/',sta1,'_',sta2,'_particle_motion.pdf'],9999,300);
            end
        end
        
    end
end