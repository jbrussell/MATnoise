% Plot the real and imaginary parts of the cross spectra to evaluate where signal 
% is best and degree of bias from inhomogeneous noise sources.
%
% https://github.com/jbrussell

clear;
setup_parameters;

%======================= PARAMETERS =======================%
comp = 'ZZ'; %'ZZ'; %'RR'; %'TT';
period_lims = [2 50];
windir = 'window3hr';
pts_smooth = 20; % just for plotting purposes
issemilogx = 0;
%==========================================================%

stalist = parameters.stalist;
nsta = parameters.nsta;
winlength = parameters.winlength;
figpath = parameters.figpath;
%fig_winlength_path = [figpath,'window',num2str(winlength),'hr/fullStack/'];
% custom directory names
    fig_winlength_path = [figpath,windir,'/fullStack/'];
dt = parameters.dt;

%------------ PATH INFORMATION -------------%
ccf_path = parameters.ccfpath;
%ccf_winlength_path = [ccf_path,'window',num2str(winlength),'hr/'];
    ccf_winlength_path = [ccf_path,windir,'/'];
ccf_singlestack_path = [ccf_winlength_path,'single/'];
ccf_daystack_path = [ccf_winlength_path,'dayStack/'];
ccf_monthstack_path = [ccf_winlength_path,'monthStack/'];
ccf_fullstack_path = [ccf_winlength_path,'fullStack/'];

ccf_stack_path = ccf_fullstack_path;

figpath = [fig_winlength_path,'bessel/'];
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
        ccf = data.coh_sum./data.coh_num;
        ccfreal = smooth(real(ccf),pts_smooth);
        ccfimag = smooth(imag(ccf),pts_smooth);
        ccfamp = smooth(abs(ccf),pts_smooth);
        
        ccf_win = data.coh_sum_win./data.coh_num;
        ccfreal_win = smooth(real(ccf_win),1);
        ccfimag_win = smooth(imag(ccf_win),1);
        ccfamp_win = smooth(abs(ccf_win),1);
        
           
        %------------------------%
        
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
        sta1sta2_dist_all(npairall) = sta1sta2_dist(nstapair); % vector containing distance between each station pair
        existpair(npairall) = {[sta1,'_',sta2]};
    
        %----------- PLOT REAL AND IMAGINARY -------------%
        f201 = figure(201); clf; set(gcf, 'Color', 'w'); box on;
        set(gcf,'position',[308   115   600   590]);
        
        subplot(3,1,1)
        dt = 1;
        T = length(ccfreal);
        faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
        ind = find(faxis>0);
        if issemilogx
            h(1) = semilogx((faxis(ind)),ccfreal(ind),'-','color',[0.5, 0.5, 0.5],'linewidth',2); hold on;
            h(2) = semilogx((faxis(ind)),ccfreal_win(ind),'-k','linewidth',2);
        else
            h(1) = plot((faxis(ind)),ccfreal(ind),'-','color',[0.5, 0.5, 0.5],'linewidth',2); hold on;
            h(2) = plot((faxis(ind)),ccfreal_win(ind),'-k','linewidth',2);
        end
        title(sprintf('REAL %s %s x-spectra %s ,distance: %f km',sta1,sta2,comp(1),sta1sta2_dist(nstapair)));
        legend(h,{'full data';'windowed'});
        xlim([1/period_lims(2) 1/period_lims(1)]);
        %xlim([0.04 0.16])
        xlabel('Frequency');
        ylim([-0.03 0.03]);
        %ylim([-0.01 0.01]);
        set(gca,'linewidth',1.5);
        
        subplot(3,1,2)
        dt = 1;
        T = length(ccfimag);
        faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
        ind = find(faxis>0);
        if issemilogx
            h(1) = semilogx((faxis(ind)),ccfimag(ind),'-','color',[0.5, 0.5, 0.5],'linewidth',2); hold on;
            h(2) = semilogx((faxis(ind)),ccfimag_win(ind),'-k','linewidth',2);
        else
            h(1) = plot((faxis(ind)),ccfimag(ind),'-','color',[0.5, 0.5, 0.5],'linewidth',2); hold on;
            h(2) = plot((faxis(ind)),ccfimag_win(ind),'-k','linewidth',2);
        end
        title(sprintf('IMAGINARY %s %s x-spectra %s ,distance: %f km',sta1,sta2,comp(1),sta1sta2_dist(nstapair)));
        xlim([1/period_lims(2) 1/period_lims(1)]);
        %xlim([0.04 0.16])
        xlabel('Frequency');
        ylim([-0.03 0.03]);
        %ylim([-0.01 0.01]);
        set(gca,'linewidth',1.5);
        
        subplot(3,1,3)
        dt = 1;
        T = length(ccfimag);
        faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
        ind = find(faxis>0);
        if issemilogx
            h(1) = semilogx((faxis(ind)),ccfamp(ind),'-','color',[0.5, 0.5, 0.5],'linewidth',2); hold on;
            h(2) = semilogx((faxis(ind)),ccfamp_win(ind),'-k','linewidth',2);
        else
            h(1) = plot((faxis(ind)),ccfamp(ind),'-','color',[0.5, 0.5, 0.5],'linewidth',2); hold on;
            h(2) = plot((faxis(ind)),ccfamp_win(ind),'-k','linewidth',2);
        end
        title(sprintf('AMP %s %s x-spectra %s ,distance: %f km',sta1,sta2,comp(1),sta1sta2_dist(nstapair)));
        xlim([1/period_lims(2) 1/period_lims(1)]);
        %xlim([0.04 0.16])
        xlabel('Frequency');
        %ylim([-0.03 0.03]);
        set(gca,'linewidth',1.5);
        
        drawnow;
%         pause;
        %print(f201,'-dpdf',[figpath,'bessel_',comp,'_',sta1,'_',sta2,'.pdf']); % Save figure
        if issemilogx
            save2pdf([figpath,'log_bessel_',comp,'_',sta1,'_',sta2,'_',num2str(pts_smooth),'pts_comparewin.pdf'],f201,1000);
        else
            save2pdf([figpath,'bessel_',comp,'_',sta1,'_',sta2,'_',num2str(pts_smooth),'pts_comparewin.pdf'],f201,1000);
        end
    end % ista2
end % ista1

%pause;