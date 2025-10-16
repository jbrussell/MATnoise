% Perform a simple FTAN analysis for estimating group velocity
% (and maximum amplitude of the Rayleigh wave for ellipticity estimates).
% The code includes some options for deciding how many mode branches to
% expect, their period ranges, and various parameters that control the
% smoothness and other characteristics of the traced group velocity curves.
% Correlations between the causal, acausal, and stacked group velocity curves
% are saved for later QC purposes. The positive, negative, and joint SNR
% values are also saved.
%
% jbrussell - 9/2025
clear all; close all;
setup_parameters;
IsFigure = 1;
IsFigure_GAUS = 1; % Plot frequency domain filtered and unfiltered
IsFigure_env = 1;

isoverwrite = 1; % overwrite results?
isoutput = 1; % save output?

%======================= PARAMETERS =======================%
%comp = {'TT', 'RR', 'ZZ'}; %'ZZ'; %'RR'; %'TT';
comp = {'ZZ'};
windir = 'window3hr';
% windir = 'window3hr_Zcorr_tiltcomp';
% Group velocity min and max
vmin = 0.5;
vmax = 6;
% Frequencies of intereset
Npers = 25; % number of periods to use
frange_fit = [1/25 1/3]; % Frequency range to estimate grv
% periods = logspace(log10(1/frange_fit(2)),log10(1/frange_fit(1)),Npers); % Period of interest for Group Velocity
periods = 1./flip(linspace(frange_fit(1),frange_fit(2),Npers)); % Period of interest for Group Velocity
%%% --- Parameters to build up gaussian filters --- %%% 
% (effects the width of the filter in the frequency domain)
% alpha = 100; % larger number = narrower bands
min_width = 0.1; %0.06; %0.18; %0.06
max_width = 0.1; %0.10; %0.30; %0.10
%%% --- Parameters forpick_ftan_ridges_robust --- %%% 
opts.Nu = 400; %           (default 400)   % # velocity samples
opts.nBranches = 2; %    (default 2)     % # dispersion branches to pick
opts.per_bounds = [6 max(periods); % vector defining period range allowed for each mode (0, 1, 2....)
                   4 8];
opts.lambda_dU = 30; %    (default 150)   % smoothness penalty on ?U (km/s)^2
opts.max_jump_U = 0.25; %   (default 0.5)   % max |?U| per period step (km/s)
opts.w_amp = 1; %        (default 1.0)   % weight on amplitude score
opts.normalizePerPeriod = true; % (true)     % robust z-score per period
opts.stopIfWeakZ = 0; % (default -Inf)  % stop if median z along path < thresh
opts.suppress_sigmaU = 0.1; % (0.12)        % Gaussian half-width (km/s) to suppress around picked ridge
opts.suppress_gain = 100; % (5)             % amount to subtract from Z along picked ridge


%==========================================================%

% dt = parameters.dt;
stalist = parameters.stalist;
nsta = parameters.nsta;
winlength = parameters.winlength;
figpath = parameters.figpath;
%fig_winlength_path = [figpath,'window',num2str(winlength),'hr/fullStack/'];
% custom directory names
fig_winlength_path = [figpath,windir,'/fullStack/'];

%------------ PATH INFORMATION -------------%
ccf_path = './ccf/';
%ccf_winlength_path = [ccf_path,'window',num2str(winlength),'hr/'];
    ccf_winlength_path = [ccf_path,windir,'/'];
ccf_singlestack_path = [ccf_winlength_path,'single/'];
ccf_daystack_path = [ccf_winlength_path,'dayStack/'];
ccf_monthstack_path = [ccf_winlength_path,'monthStack/'];
ccf_fullstack_path = [ccf_winlength_path,'fullStack/'];

ccf_stack_path = ccf_fullstack_path;


% output path
ftan_out_path = ['./ftan/',windir,'/fullStack/ftan',comp{1},'/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
if ~exist(ftan_out_path)
    mkdir(ftan_out_path)
end

% figure output path
ftan_fig_path = ['./figs/',windir,'/fullStack/ftan/',comp{1},'_',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_br',num2str(opts.nBranches),'/'];
if ~exist(ftan_fig_path)
    mkdir(ftan_fig_path);
end

ccf_path = [ccf_stack_path,'ccf',comp{1},'/',];
npairall = 0;
%------------ LOAD DATA AND PLOT IN TIME DOMAIN -------------%
for ista1=1:nsta % loop over all stations
    sta1=char(stalist(ista1,:));
    sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
    
    nstapair = 0;
    ccf_filt_env_sum = 0;
    for ista2 = 1: nsta % loop over station pairs
        sta2 = char(stalist(ista2,:));
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        
        if ~isoverwrite
            if exist(sprintf('%s/%s_%s_ftan.mat',ftan_out_path,sta1,sta2))
                disp(['Already processed ',sta1,'-',sta2,' ... skipping']);
                continue;
            end
        end
                
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        
        if ~exist(filename,'file') % check that ccf file exists
            disp(['not exist ',filename])
            continue;
        end
        nstapair = nstapair + 1;
        sta2_list{nstapair} = sta2;
        
        %----------- LOAD DATA -------------%
        data = load(filename);
        dt = data.stapairsinfo.dt;
        ccf = data.coh_sum./data.coh_num;
        ccf_orig = ccf; % save original unfiltered CCF
        
        % Distance between sta1 and sta2
        r = distance(data.stapairsinfo.lats(1),data.stapairsinfo.lons(1),data.stapairsinfo.lats(2),data.stapairsinfo.lons(2),referenceEllipsoid('GRS80'))/1000;

        %% Apply Gaussian filters to data in frequency domain
        
        [ccf_gaus,faxis, gaus_filters] = gaus_filt_nbands(ccf,periods,dt,min_width,max_width);
%         [ccf_gaus,faxis, gaus_filters] = gaus_filt_nbands_levshin(ccf,periods,dt,alpha);
        
        % Calculate SNR for positive, negative, and both
        isfigure_snr = 0;
        snr = zeros(size(periods));
        snr_pos = zeros(size(periods));
        snr_neg = zeros(size(periods));
        for iper = 1:size(ccf_gaus,2)
            [snr(iper), signal_ind, snr_pos(iper), snr_neg(iper)] = calc_SNR(ccf_gaus(:,iper),vmin,vmax,r,dt,isfigure_snr);
%             pause;
        end
    
        % Plot the gaussian filters and an example of a filtered applied in the frequency domain
        if IsFigure_GAUS

            figure(1); clf;
            subplot(2,1,1);
            [xi yi] = ndgrid(faxis,1./periods);
            surface(xi,yi,gaus_filters);
            shading flat;
            xlabel('Frequency axis (Hz)');
            ylabel('Narrow band filter freq (Hz)');
        % 	xlim([0 max(centf)*1.5])
            title(['Gaussian filters ',sta1,'-',sta2]);
            set(gca,'linewidth',1.5,'fontsize',15);
            axis tight;
            
            subplot(2,1,2); box on; hold on;
            iper = 10;
            plot(faxis,abs(ccf),'-k','linewidth',3); hold on;
            plot(faxis,abs(ccf_gaus(:,iper)),'r','linewidth',2);
            xlabel('Frequency (Hz)');
            ylabel('Spectral Power');
            title(['Example narrow band filter: ',num2str(1./periods(iper)),' Hz']);
            legend({'Raw','Filtered'},'location','northeast');
            set(gca,'linewidth',1.5,'fontsize',15);
            axis tight
            
        end
        
        %----------- Frequency ==> Time domain -------------%
        N = size(ccf_gaus,1);
        ccf_ifft = real(ifft(ccf_gaus,N)); % inverse FFT to get time domain
        ccf_ifft = fftshift(ccf_ifft,1); % rearrange values as [-lag lag]
        %rearrange and keep values corresponding to lags: -(len-1):+(len-1)
        ccf_orig_ifft = real(ifft(ccf_orig,N)); % inverse FFT to get time domain
        ccf_orig_ifft = fftshift(ccf_orig_ifft); % rearrange values as [-lag lag]
        
        % Plot some of the filtered CCFs just to check
        if IsFigure_env
            f2 = figure(2); clf;
            set(gcf,'position',[ 60    79   743   939]);
            box on; hold on;
            time = ([0:N-1]-floor(N/2))*dt;  % build lagtime vector for plotting
            time = [time(time<0), time(time>=0)];
            N = size(ccf_gaus,1);
            plot(time,ccf_orig_ifft./max(abs(ccf_orig_ifft)),'-k','linewidth',1.5); hold on;
            text(500*1.05,1,'Unfiltered','fontsize',12,'fontweight','bold');
            xlim([-500 500]);
            clrs = jet(length(periods));
            icnt = 0;
            for iper = 1:floor(length(periods)/10):length(periods)
                icnt = icnt+1;
                plot(time,ccf_ifft(:,iper)./max(abs(ccf_ifft(:,iper)))-2*icnt,'-','linewidth',1.5,'color',clrs(iper,:));
                text(500*1.05,-2*icnt,[num2str(periods(iper)),' s'],'fontsize',12,'fontweight','bold');
                text(-500*1.05,-2*icnt,['SNR=',num2str(snr(iper))],'fontsize',12,'fontweight','bold','HorizontalAlignment','right');
            end
            set(gca,'fontsize',15,'linewidth',1.5);
            title(['Gaussian filtered CCFs: ',sta1,'-',sta2]);
            xlabel('Lag Time (s)');
%             ylabel('Normalized & Shifted Amplitudes');
        end
        
        % Get time and group velocity axes
        time = ([0:N-1]-floor(N/2))*dt;  % build lagtime vector for plotting
        time = [time(time<0), time(time>=0)];
        grv_axis = abs(r./time(time>=0)); % index only positive time
        
        % Calculate analytic envelope of waveform
        ccf_ifft_env = abs(hilbert(ccf_ifft)); % hilbert transform to get analytic envelope for each narrow band filter
%         ccf_ifft_env = log(ccf_ifft_env);

        %% Pick group velocity times

        % Pick group times by ridge tracking with constraints (causal +times)
        ind_tpos = find(time>=0);
        ccf_ifft_env_pos = ccf_ifft_env(ind_tpos,:);
        time_pos = time(ind_tpos);
        [tg_pos, grv_pos, amp_pos, idx_tg_pos, meta_pos] = pick_ftan_ridges_robust(ccf_ifft_env_pos, time_pos, periods, r, vmin, vmax,opts);
        
        % Pick group times by ridge tracking with constraints (acausal -times)
        ind_tneg = find(time<=0);
        ccf_ifft_env_neg = flip(ccf_ifft_env(ind_tneg,:));
        time_neg = flip(abs(time(ind_tneg)));
        [tg_neg, grv_neg, amp_neg, idx_tg_neg, meta_neg] = pick_ftan_ridges_robust(ccf_ifft_env_neg, time_pos, periods, r, vmin, vmax, opts);
        
        % Pick group times by ridge tracking with constraints (Stack + and -)
        ind_tneg = find(time<=0);
        ccf_ifft_env_stack = 0.5*(ccf_ifft_env_pos+ccf_ifft_env_neg);
        [tg_stack, grv_stack, amp_stack, idx_tg_stack, meta_stack] = pick_ftan_ridges_robust(ccf_ifft_env_stack, time_pos, periods, r, vmin, vmax, opts);
        
        %% Calculate correlation coefficients for pos, neg, and stack
        
        if IsFigure
            figure(4); clf;
            set(gcf,'position',[527          40        1065         420],'color','w');
            corr_mat = {};
            corr_avg = [];
            for ibr = 1:opts.nBranches 
                [R,P,RLO,RUP] = corrcoef([grv_pos(:,ibr),grv_neg(:,ibr),grv_stack(:,ibr)],'rows', 'complete');
                corr_mat{ibr} = R;
                corr_avg(ibr) = mean(unique(R(R<1)));
                
                subplot(1,opts.nBranches,ibr);
                colormap(jet);
                imagesc(R);
                cb = colorbar;
                ylabel(cb,'Correlation Coeff.');
                set(cb,'linewidth',1.5);
                caxis([0.7 1]);
                set(gca,'fontsize',15,'linewidth',1.5,'YDir','reverse');
                yticks([1 2 3]);
                xticks([1 2 3]);
                yticklabels({'+Lag';'-Lag';'Stack'})
                xticklabels({'+Lag';'-Lag';'Stack'})
                ylim([0.5 3.5]);
                xlim([0.5 3.5]);
                axis square;
                title(['Branch: ',num2str(ibr),'   R_{av}=',num2str(corr_avg(ibr))]);
            end
            sgtitle('Dispersion Curve Correlations','fontsize',16,'fontweight','bold');
            if isoutput
                save2pdf([ftan_fig_path,'/',sta1,'_',sta2,'_correlations.pdf'],4,300);
            end
        end
        
        %% Plot group velocity measurements
        if IsFigure
            f3 = figure(3); clf; 
            set(gcf, 'Color', 'w','position',[681         124        1055         420*2]);
            cmap = parula;
            colormap(cmap)
            [PERIODS, GRV] = meshgrid(periods,grv_axis);
            
            clims = [-3 0];
            
            % CAUSAL (+LAG)
            subplot(2,2,2); box on; hold on;
            set(gca,'color',cmap(1,:))
            levels = linspace(clims(1),clims(2),25);
            contourf(PERIODS,GRV,log(ccf_ifft_env_pos./max(ccf_ifft_env_pos(:))),levels,'LineStyle','none'); shading flat;
            plot(periods,grv_pos,'o-','color','r','linewidth',1.5);
            cb = colorbar;
            ylabel(cb,'Log(Amp.)');
            set(cb,'linewidth',1.5);
    %         caxis([0 1]);
            caxis(clims);
            axis tight;
            xlabel('Period (s)','fontsize',15);
            ylabel('Group Velocity (km/s)','fontsize',15);
            set(gca,'fontsize',15);
            ylim([vmin,vmax]);
            xlim([periods(1) periods(end)]);
            title(['Causal ',sta1,'-',sta2,': ',num2str(r),'km']);
            
            % ACAUSAL (-LAG)
            subplot(2,2,1); box on; hold on;
            set(gca,'color',cmap(1,:))
            contourf(PERIODS,GRV,log(ccf_ifft_env_neg./max(ccf_ifft_env_neg(:))),levels,'LineStyle','none'); shading flat;
            plot(periods,grv_neg,'o-','color','r','linewidth',1.5);
            cb = colorbar;
            ylabel(cb,'Log(Amp.)');
            set(cb,'linewidth',1.5);
    %         caxis([0 1]);
            caxis(clims);
            axis tight;
            xlabel('Period (s)','fontsize',15);
            ylabel('Group Velocity (km/s)','fontsize',15);
            set(gca,'fontsize',15);
            ylim([vmin,vmax]);
            xlim([periods(1) periods(end)]);
            title(['Acausal ',sta1,'-',sta2,': ',num2str(r),'km']);
            
            % STACKED (POS + NEG)
            subplot(2,2,3); box on; hold on;
            set(gca,'color',cmap(1,:))
            contourf(PERIODS,GRV,log(ccf_ifft_env_stack./max(ccf_ifft_env_stack(:))),levels,'LineStyle','none'); shading flat;
            plot(periods,grv_stack,'o-','color','r','linewidth',1.5);
            cb = colorbar;
            ylabel(cb,'Log(Amp.)');
            set(cb,'linewidth',1.5);
    %         caxis([0 1]);
            caxis(clims);
            axis tight;
            xlabel('Period (s)','fontsize',15);
            ylabel('Group Velocity (km/s)','fontsize',15);
            set(gca,'fontsize',15);
            ylim([vmin,vmax]);
            xlim([periods(1) periods(end)]);
            title(['Stacked ',sta1,'-',sta2,': ',num2str(r),'km']);
            
            if isoutput
                save2pdf([ftan_fig_path,'/',sta1,'_',sta2,'_grv_panels.pdf'],f3,300);
            end
        end
        
        
        %% Save outputs
        
        if isoutput
            % Package output
            ftan.periods = periods;
            ftan.grv_axis = grv_axis;
            ftan.time = time;
            ftan.snr = snr;
            for ibr = 1:opts.nBranches
                ftan.pos.grv = grv_pos;
                ftan.pos.tg = tg_pos;
                ftan.pos.amp = amp_pos;
                ftan.pos.idx_tg = idx_tg_pos;
                ftan.pos.ccf_envs = ccf_ifft_env_pos;
                ftan.pos.snr = snr_pos;
                ftan.pos.meta = meta_pos;

                ftan.neg.grv = grv_neg;
                ftan.neg.tg = tg_neg;
                ftan.neg.amp = amp_neg;
                ftan.neg.idx_tg = idx_tg_neg;
                ftan.neg.ccf_envs = ccf_ifft_env_neg;
                ftan.neg.snr = snr_neg;
                ftan.neg.meta = meta_neg;

                ftan.stack.grv = grv_stack;
                ftan.stack.tg = tg_stack;
                ftan.stack.amp = amp_stack;
                ftan.stack.idx_tg = idx_tg_stack;
                ftan.stack.ccf_envs = ccf_ifft_env_stack;
                ftan.stack.meta = meta_stack;
            end
            ftan.corr_mat = corr_mat;
            ftan.corr_avg = corr_avg;
            ftan.vmin = vmin;
            ftan.vmax = vmax;
            ftan.min_width = min_width;
            ftan.max_width = max_width;
            ftan.r = r;
            ftan.sta1 = sta1;
            ftan.sta2 = sta2;
            ftan.stapairsinfo = data.stapairsinfo;
            
            save(sprintf('%s/%s_%s_ftan.mat',ftan_out_path,sta1,sta2),'ftan');
        end

        
    end % ista2
end % ista1
