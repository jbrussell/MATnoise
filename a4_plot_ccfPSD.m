% Plot the power spectral density to identify frequency content of signal
%
% https://github.com/jbrussell

clear;
setup_parameters;

%======================= PARAMETERS =======================%
comp = 'ZZ'; %'ZZ'; %'RR'; %'TT';
windir = 'window3hr'; 
xlim_per = [1/10 100]; % Period bounds for plotting

%==========================================================%

stalist = parameters.stalist;
nsta = parameters.nsta;
winlength = parameters.winlength;
figpath = parameters.figpath;
%fig_winlength_path = [figpath,'window',num2str(winlength),'hr/fullStack/'];
% custom directory names
    fig_winlength_path = [figpath,windir,'/fullStack/'];

%------------ PATH INFORMATION -------------%
ccf_path = parameters.ccfpath;
%ccf_winlength_path = [ccf_path,'window',num2str(winlength),'hr/'];
    ccf_winlength_path = [ccf_path,windir,'/'];
ccf_singlestack_path = [ccf_winlength_path,'single/'];
ccf_daystack_path = [ccf_winlength_path,'dayStack/'];
ccf_monthstack_path = [ccf_winlength_path,'monthStack/'];
ccf_fullstack_path = [ccf_winlength_path,'fullStack/'];

ccf_stack_path = ccf_fullstack_path;

figpath = [fig_winlength_path,'psd/'];
% create figure directory
if ~exist(fig_winlength_path)
    mkdir(fig_winlength_path)
end
if ~exist(figpath)
    mkdir(figpath)
end

ccf_path = [ccf_stack_path,'ccf',comp,'/',];
npairall = 0;
sum_all_psd = 0;
sum_all_psd_log = 0;
%------------ LOAD DATA AND CALCULATE AND PLOT PSD -------------%
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
        dt = data.stapairsinfo.dt;
        ccf = data.coh_sum_win./data.coh_num;
        if ~isempty(ccf(isnan(ccf)))
            display(['skipping ',sta1,'-',sta2]);
            continue
        end
        
        %----------- CALCULATE PSD -------------%
        N = length(ccf);
        tmax = N*dt;
        X = real(ifft(ccf)); % time domain
        Fs = 1/dt;
        WINDOW = []; %[round(length(ccf)/10)];
        NOVERLAP = [];
        NFFT = [];
        [ccf_psd{nstapair},F] = pwelch(X,WINDOW,NOVERLAP,NFFT,Fs,'onesided');
        ccf_psd_log{nstapair} = 10*log10(ccf_psd{nstapair});
        %ccf_psd{nstapair} = (2/tmax)*abs(ccf*dt).^2;
        %pts_smooth = 10;
        %ccf_psd_log{nstapair} = smooth(10*log10(ccf_psd{nstapair}),pts_smooth);
%         ccf_psd{nstapair} = (1/N/Fs) * abs(ccf(1:floor(N/2+1))).^2 * 2;
%         ccf_psd_log{nstapair} = smooth(10*log10(ccf_psd{nstapair}),100);

        
        % Distance between sta1 and sta2
        sta1sta2_dist(nstapair) = deg2km(distance(data.stapairsinfo.lats(1),data.stapairsinfo.lons(1),data.stapairsinfo.lats(2),data.stapairsinfo.lons(2)));
        
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
        ccf_psd_all{npairall} = ccf_psd{nstapair} ; % cell containing all psd
        ccf_psd_log_all{npairall} = ccf_psd_log{nstapair};
        sta1sta2_dist_all(npairall) = sta1sta2_dist(nstapair); % vector containing distance between each station pair
        existpair(npairall) = {[sta1,'_',sta2]};
        sum_all_psd = sum_all_psd + ccf_psd{nstapair};
        sum_all_psd_log = sum_all_psd_log + ccf_psd_log{nstapair};
    end % ista2
    
        %----------- CALCULATE MEAN STATION PSD  -------------%
        sum_sta_psd = 0;
        sum_sta_psd_log = 0;
        nsta2 = length(ccf_psd);
        for ista2 = 1:nsta2
            if isempty(ccf_psd{ista2})
                continue
            end
            sum_sta_psd = sum_sta_psd + ccf_psd{ista2};
            sum_sta_psd_log = sum_sta_psd_log + ccf_psd_log{ista2};
        end
        ccf_psd_staMean{ista1} = sum_sta_psd/nsta2;
        ccf_psd_log_staMean{ista1} = sum_sta_psd_log/nsta2;
    
        %----------- PLOT STATION PSDs -------------%
        f201 = figure(201); clf; hold on; set(gcf, 'Color', 'w'); box on;
        T = length(ccf);
        faxis = F;
        ind = find(faxis>=0);
        clr = lines(nstapair);
        
        subplot(2,1,1); % plot wide period band
        for istapair = 1: nstapair % loop over station pairs
            if isempty(ccf_psd_log{istapair})
                continue
            end
            semilogx((1./faxis(ind)),ccf_psd_log{istapair}(ind),'-','color',[.5 .5 .5]); hold on;
            %pause
        end
        h1 = semilogx((1./faxis(ind)),ccf_psd_log_staMean{ista1}(ind),'-k','linewidth',3);
        xlim(xlim_per);
        xlabel('Period (s)','fontsize',18);
        ylabel('Power','fontsize',18);
        title(['reference station:',sta1,' ',comp(1)],'fontsize',18,'fontweight','bold');
        legend(h1,{'mean'},'fontsize',12);
        set(gca,'fontsize',15);
        
        subplot(2,1,2); % plot short periods only
        for istapair = 1: nstapair % loop over station pairs
            if isempty(ccf_psd_log{istapair})
                continue
            end
            plot((1./faxis(ind)),ccf_psd_log{istapair}(ind),'-','color',[.5 .5 .5]); hold on;
            %pause
        end
        h1 = plot((1./faxis(ind)),ccf_psd_log_staMean{ista1}(ind),'-k','linewidth',3);
        xlim(xlim_per);
        xlabel('Period (s)','fontsize',18);
        ylabel('Power','fontsize',18);
        title(['reference station:',sta1,' ',comp(1)],'fontsize',18,'fontweight','bold');
        legend(h1,{'mean'},'fontsize',12);
        set(gca,'fontsize',15);

end % ista1

%----------- CALCULATE MEAN NETWORK PSD -------------%
ccf_psd_allMean = sum_all_psd/npairall;
ccf_psd_log_allMean = sum_all_psd_log/npairall;

%----------- PLOT NETWORK AVERAGE PSD -------------%
f203 = figure(203); clf; hold on; set(gcf, 'Color', 'w'); box on;
clr = lines(npairall);

subplot(2,1,1); % Plot wide period band
for istapair = 1: npairall
    if isempty(ccf_psd_log_all{istapair})
        continue
    end
    semilogx((1./faxis(ind)),ccf_psd_log_all{istapair}(ind),'-','color',[.5 .5 .5]); hold on;
end
h1 = semilogx((1./faxis(ind)),ccf_psd_log_allMean(ind),'-k','linewidth',3);
axis tight;
xlim(xlim_per);
ylim([-140 -80]); 
xlabel('Period (s)','fontsize',18);
ylabel('Power','fontsize',18);
title(['All stations ',comp(1)],'fontsize',18,'fontweight','bold');
set(gca,'fontsize',15,'linewidth',1.5);
set(gca,'Xtick',[1 2 3 4 5 6 7 8 10 20 30 40 50 70 100]);

%pause;

% print(f203,'-dpdf',[figpath,'psd_ccfwin_allsta',comp,'_log_pwelch.pdf']); % Save figure
save2pdf([figpath,'psd_ccfwin_allsta',comp,'_log_pwelch_ORALS.pdf'],f203,1000);
% export_fig([figpath,'psd_ccfwin_allsta',comp,'_log_pwelch_ORALS.pdf'],'-pdf','-q100','-p0.02','-painters',f203)
