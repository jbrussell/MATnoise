% Extract phase velocity dispersion between station pairs by fitting J0 bessel 
% function to real(ccf)
% Uses cross spectral fitting technique of Menke & Jin (2015) BSSA 
% DOI:10.1785/0120140245
%
% Define own starting phase velocity dispersion c manually or using 
% functions/calc_Rayleigh_disp for a simple layered model (does not work for 
% models with a water column).
%
%
% https://github.com/jbrussell
clear
close all;

global tN
global waxis
global twloc
global weight
setup_parameters;

%======================= PARAMETERS =======================%

% RAYLEIGH FUND MODE
comp = {'ZZ'}; %'RR'; 'ZZ'; 'TT'
windir = 'window3hr';
xspdir = 'ZZ_0S_LRT'; % output directory of phase velocities
N_wl = 1; % Number of wavelengths required
Ninterp = 40; % [] or Number of points to interpolate to;
% Use picks from Linear Radon Transform? (./mat-LRTdisp/)
is_LRT_picks = 1; % Use picks from Radon Transform to determine starting dispersion model and frequencies
LRT_method = 'CGG_weight';
mode_br = 0; % desired mode branch [0=fund.]
frange_LRT = [1/40 1/3]; % Frequency range of LRT panel for reading in picks
frange_fit = [1/25 1/3]; % Frequency range to fit over! Can be more restrictive than where picks were made
damp = [1; 1; 1]; % [fit, smoothness, slope]
is_normbessel = 0; % normalize bessel function by analytic envelope?

if ~is_LRT_picks
    frange_fit = [1/40 1/10]; % frequency range over which to fit bessel function
%     xlims = [1/70 1/9];
    Npers = 18; % Number of periods
    t_vec_all = 1./flip(linspace(frange_fit(1) , frange_fit(2) ,Npers)); % periods at which to extract phase velocity
end

xlims = [min(frange_fit)*0.9 max(frange_fit)*1.1]; % limits for plotting

is_resume = 0; % Resume from last processed file (1) or overwrite (0)
iswin = 1; % Use the time-domain windowed ccfs?
npts_smooth = 1; % Smoothing of bessel function: 1 = no smoothing

isoutput = 1; % Save *.mat file with results?
nearstadist = 0;
IsFigure = 1;
isfigure2 = 0;
isfigure_snr = 1;

%% Make initial guess at phase velocity dispersion model
if is_LRT_picks
    % Read picks from Linear Radon Transform (see ./mat-LRTdisp/)
    ccfstr = strsplit(parameters.ccfpath,'/');
    ccfstr = ccfstr{end-1};
    in_LRTpicks = [parameters.path_LRT_picks,ccfstr,'/',windir,'/',num2str(1/frange_LRT(2)),'_',num2str(1/frange_LRT(1)),'s/LRTpicks_',LRT_method,'_',comp{1},'.mat'];
    temp = load(in_LRTpicks);
    picks_LRT = temp.picks_LRT;
    
    imode = mode_br+1;
     % Interpolate to number of desired points
    if ~isempty(Ninterp)
        t_vec_int = 1./(linspace( 1./min(picks_LRT(imode).per), 1./max(picks_LRT(imode).per), Ninterp ));
        picks_LRT(imode).phv = interp1( picks_LRT(imode).per, picks_LRT(imode).phv, t_vec_int );
        picks_LRT(imode).phv_std = interp1( picks_LRT(imode).per, picks_LRT(imode).phv_std, t_vec_int );
        picks_LRT(imode).per = t_vec_int;
    end
    
    % Avoid overlapping periods from other mode_br branches
    if mode_br == 0
        if imode < length(picks_LRT)
            I_good = picks_LRT(imode).per > picks_LRT(imode+1).per(end);
        else
            I_good = true(size(picks_LRT(imode).per));
        end
    elseif mode_br > 0
        if imode-1 ~= 0 && imode < length(picks_LRT)
            I_good = (picks_LRT(imode).per < picks_LRT(imode-1).per(1)) & (picks_LRT(imode).per > picks_LRT(imode+1).per(end));
        end
        if imode-1 ~= 0
            I_good = picks_LRT(imode).per < picks_LRT(imode-1).per(1);
        end
        if imode < length(picks_LRT)
            I_good = picks_LRT(imode).per > picks_LRT(imode+1).per(end);
        end
    end
    
    % Fit only periods within frange_fit
    I_fit = picks_LRT(imode).per>=1/frange_fit(2) & picks_LRT(imode).per<=1/frange_fit(1);
    I_good = logical(I_fit .* I_good);
    
    c_all = picks_LRT(imode).phv(I_good);
    t_vec_all = picks_LRT(imode).per(I_good);
    c_all_std = picks_LRT(imode).phv_std(I_good);
    c_start = c_all;
else
    % Read from MINEOS .q file (https://github.com/jbrussell/MINEOS_synthetics)
    qfile = ['./qfiles/Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_ORCAiso_INV.s0to200.q'];
    if exist('c','var') == 0 % check if phase velocities exist, if not read them in
        [~,~,c_all] = readMINEOS_qfile2(qfile,t_vec_all,mode_br);
    end
    c_start = c_all;
    c_all_std = zeros(size(c_all));
end

%%
%==========================================================%

% LIMITS
if comp{1}(1) == 'R'
    ylims = [3.2 4.5];
elseif comp{1}(1) == 'Z' || comp{1}(1) == 'P'
%     ylims = [1.5 5.0];
%     ylims = [2 4.5];
    ylims = [1.5 4.5];
elseif comp{1}(1) == 'T'
    ylims = [3.5 4.8];
end

% LOAD DATA TO SEE HOW MANY POINTS
%%% --- Load in the ccf --- %%%
        %ccf_path = ['./ccf/',windir,'/fullStack/ccf',comp{1},'/'];
        ccf_path = [parameters.ccfpath,windir,'/fullStack/ccf',comp{1},'/'];
        stalist = parameters.stalist;
        sta1=char(stalist(1,:));
        sta2=char(stalist(2,:));
        sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        if ~exist(filename,'file')
            disp(['not exist ',filename])
%             continue;
        end
        data1 = load(filename);
        npts = length(data1.coh_sum_win);


% input path
%ccf_path = ['./ccf/',windir,'/fullStack/ccf',comp{1},'/'];
ccf_path = [parameters.ccfpath,windir,'/fullStack/ccf',comp{1},'/'];

% output path
XSP_path = ['./Xsp/',windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_',num2str(N_wl),'wl_',xspdir,'/'];

if ~exist(XSP_path)
    if ~exist('./Xsp/')
        mkdir('./Xsp/');
    end
    if ~exist(['./Xsp/',windir,'/'])
        mkdir(['./Xsp/',windir,'/']);
    end
    if ~exist(['./Xsp/',windir,'/fullStack/'])
        mkdir(['./Xsp/',windir,'/fullStack/']);
    end
    if ~exist(['./Xsp/',windir,'/fullStack/Xsp',comp{1},'/'])
        mkdir(['./Xsp/',windir,'/fullStack/Xsp',comp{1},'/']);
    end
    mkdir(XSP_path)
end

% figure output path
if iswin
    XSP_fig_path = ['./figs/',windir,'/fullStack/Xsp/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_',num2str(N_wl),'wl_',xspdir,'/TEI19/'];
else
    XSP_fig_path = ['./figs/',windir,'/fullStack/Xsp/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_',num2str(N_wl),'wl_',xspdir,'/TEI19_nowin/'];
end

if ~exist(XSP_fig_path)
    mkdir(XSP_fig_path);
end




warning off; %#ok<WNOFF>

stalist = parameters.stalist;
nsta=parameters.nsta; % number of target stations to calculate for

%%% --- Loop through station 1 --- %%%
for ista1=1:nsta
    
    sta1=char(stalist(ista1,:));
    sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
    
    %%% --- Loop through station 2 --- %%%
    for ista2 = 1: nsta % length(v_sta)
        sta2 = char(stalist(ista2,:));
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        
        % Check to see if we have already done this
        if is_resume && exist([XSP_path,sta1,'_',sta2,'_xsp.mat'])
            disp('Already fit this one!')
            continue
        end
        clear data1 xcorf1 xsp1 filename
        
        %%% --- Load in the ccf --- %%%
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        
        if ~exist(filename,'file')
            disp(['not exist ',filename])
            continue;
        end
        
        
        data1 = load(filename);
        r1 = distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2),referenceEllipsoid('GRS80'))/1000;
        dt = data1.stapairsinfo.dt;
        groupv_max = data1.max_grv;
        groupv_min = data1.min_grv;
        
        if r1 < nearstadist
            continue;
        end
        
        % Index wavelength criterion
        I_wl = r1 ./ (t_vec_all .* c_all) > N_wl;
        if sum(I_wl) <= 1
            I_wl(1) = 1;
            I_wl(2) = 1;
        end
        c = c_all(I_wl);
        t_vec = t_vec_all(I_wl);
        c_std = c_all_std(I_wl);
        
        tN = length(t_vec);
        wholesec = npts;
        wvec1 = (2*pi)./t_vec;
        wvec1 = wvec1';
        
        % Get your axis correct
        twloc=1./t_vec;
        twloc = twloc*2*pi;
%         waxis = (frange_fit(1):1/wholesec:frange_fit(2))*2*pi;
        waxis = (1/max(t_vec):1/wholesec:1/min(t_vec))*2*pi;
        
        
        %%% - Get the normalized ccf - %%%
        if iswin
            xcorf1 = data1.coh_sum_win./data1.coh_num;
        else
            xcorf1 = data1.coh_sum./data1.coh_num;
        end
        dumnan = find(isnan(xcorf1)==1);
        
        if length(dumnan) > 10
            disp([sta1,' and ',sta2,'is NaN! Moving on']);
            continue
        end
        
        %N = 10000;
        N = length(xcorf1);
        if length(xcorf1) < N
            disp('Dataset is too short! Moving on')
            continue
        end
        
        xcorf = xcorf1;
        xcorf1 = real(xcorf1(1:N));
        xcorf1(1) = 0;
        
        if isfigure2 
            figure(1)
            T = length(xcorf1);
            temp_faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(temp_faxis>0);
            subplot(2,1,1)
            %plot(temp_faxis(ind),smooth(real(xcorf1(ind)),100));
            plot(flip(temp_faxis(ind),smooth(real(xcorf1(ind)),50)));
            xlim([frange_fit(1) frange_fit(2)])
            hold on
            subplot(2,1,2)
            %plot(temp_faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),100),'-r')
            plot(temp_faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),50),'-r');
            xlim([frange_fit(1) frange_fit(2)])
            
        end

        %%% - Convert xcorf into spherical frequency - %%%
        faxis = [0:(N-mod(N-1,2))/2 , -(N-mod(N,2))/2:-1]/dt/N;
        xsp1 = interp1(faxis*2*pi,xcorf1,waxis);

        %xsp1 = smooth(xsp1,50);
        xsp1 = smooth(xsp1,npts_smooth);

        tw1 = ones(1,tN)*r1./c;
        
        %%% - Invert for the bessel function 2x - %%%
        options = optimoptions(@lsqnonlin,'TolFun',1e-12,'MaxIter',1500,'MaxFunEvals',1500);
        weight  = 1./waxis;
        tw2 = lsqnonlin(@(x) besselerr(x,[xsp1],damp,is_normbessel),[tw1],[tw1]*0.8,[tw1]*1.2,options);
%         tw2 = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw1],[],[],options);
        
        weight(:) = 1;
        [tw,~,res,~,~,~,J] = lsqnonlin(@(x) besselerr(x,[xsp1],damp,is_normbessel),[tw2],[tw2]*0.8,[tw2]*1.2,options);
%         tw = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw2],[tw2]*0.8,[tw2]*1.2,options);
%         tw = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw2],[],[],options);
        
        % ESTIMATE ERROR BARS ON MODEL! :: JBR - 2/2020
        % Calculate data variance from residual following Menke QMDA book eq. (4.31)
        % s_d^2 = E / (N-M)
        sigma_d2 = res'*res / (length(res)-length(tw));
        Cov_m = inv(J'*J)*sigma_d2;
        sigma_m_tw = diag(Cov_m).^(1/2);
        % Propogate model error to phase velocity
        % c = r1./tw;  therefore   dc = |r*t^(-2) * dt|
        sigma_m_c = abs(r1.*tw'.^(-2).*sigma_m_tw);
        
        %%% - Set up the variable structure - %%%
        xspinfo.sta1 = sta1;
        xspinfo.sta2 = sta2;
        xspinfo.lat1 = data1.stapairsinfo.lats(1);
        xspinfo.lon1 = data1.stapairsinfo.lons(1);
        xspinfo.lat2 = data1.stapairsinfo.lats(2);
        xspinfo.lon2 = data1.stapairsinfo.lons(2);
        
        xspinfo.r = r1;
        xspinfo.tw = tw;
        xspinfo.xsp = xsp1;
        xspinfo.xsp_norm = xsp1./abs(hilbert(xsp1));
        xspinfo.coherenum = data1.coh_num;
        err = besselerr(tw,xsp1,damp,is_normbessel);
        err = err(1:length(waxis));
        
        if is_normbessel
            xspinfo.sumerr = sum(err.^2)./sum((xspinfo.xsp_norm./weight(:)).^2);
        else
            xspinfo.sumerr = sum(err.^2)./sum((xsp1./weight(:)).^2);
        end
        xspinfo.err = err./weight(:);
        xspinfo.tw1 = tw1;
        xspinfo.twloc = twloc;
        xspinfo.c = r1./tw;
        xspinfo.c_std = sigma_m_c;
        xspinfo.per = 1./(twloc/2/pi);
        xspinfo.c_start = c_start;
        xspinfo.c_std_start = c_all_std;
        xspinfo.per_start = t_vec_all;
        xspinfo.isgood_wl = I_wl;
        
        data = r1./tw;
        

        %% %%% Calculate SNR %%%
        xcorf1 = data1.coh_sum./data1.coh_num;
        xcorf1_filtered = tukey_filt( xcorf1,[min(t_vec) max(t_vec)],dt,0.25 );
        [snr, signal_ind] = calc_SNR(xcorf1_filtered,groupv_min,groupv_max,r1,dt,isfigure_snr);
        %%

        xspinfo.filename = filename;
        xspinfo.snr = snr;
        
        % Calculate the predicted bessel function from the initial model
        disp([filename,' fitted'])
        if IsFigure
            if 0 % plot initial bessel
                tw_init = interp1(twloc,tw1(1:tN),waxis,'linear');
                x_init = waxis.*tw_init;
                A = 1;
                binit = besselj(0,x_init)*A;
                binit = binit./mean(abs(binit)).*mean([abs(xsp1)]);
                plot(waxis/2/pi,binit,'-k','linewidth',2); hold on;
            end
            
            f3 = figure(3); clf; hold on; 
            set(gcf,'color','w','Position',[289     1   517   704]);

            ax1 = subplot(3,1,1);
            plot_SNR(xcorf1_filtered,groupv_min,groupv_max,r1,dt,ax1);
            set(ax1,'box','off');

            % REAL PART (J0)
            subplot(3,1,2); box off;
            tww = interp1(twloc,tw(1:tN),waxis,'linear');
            x = waxis.*tww;
            A = 1;
            b = besselj(0,x)*A;
            b = b./mean(abs(b)).*mean([abs(xsp1)]);           
            T = length(data1.coh_sum);
            faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(faxis>0);
            plot(faxis(ind),smooth(real(data1.coh_sum_win(ind)/data1.coh_num),npts_smooth),'-k','linewidth',3); hold on;
            if ~iswin
                plot(waxis/2/pi,xsp1,'-b','linewidth',1);
            end
            plot(waxis/2/pi,b,'-r','linewidth',2); hold on; 
            xlim(xlims);
            xlims1 = get(gca,'XLim');
            ylabel('J_{0}','fontsize',16);
            ax1 = get(gca);
            dx = 1;
            dy = 0.95;
            text(gca,diff(ax1.XLim)*dx+ax1.XLim(1),diff(ax1.YLim)*dy+ax1.YLim(1),...
                sprintf('%.4f',xspinfo.sumerr),'color',[0 0 0],'fontsize',14);
            set(gca,'fontsize',16,'linewidth',1.5);
            box off;
            title(['Distance : ',num2str(r1),'km'],'fontsize',16);
                
            hold on
            subplot(3,1,3);
            errorbar(twloc/2/pi,r1./tw1,c_std,'o-','color',[0.5 0.5 0.5],'linewidth',2);hold on;
            errorbar(twloc/2/pi,r1./tw,sigma_m_c*2,'o-','color',[1 0 0 ],'linewidth',2);
%             errorbar(twloc/2/pi,r1./tw,xspinfo.err,'ro-','linewidth',2);
            title([sta1,'-',sta2],'fontsize',16)
            xlabel('Frequency (Hz)','fontsize',16);
            ylabel('Phase Velocity (km/s)','fontsize',16);
            set(gca,'fontsize',16,'linewidth',1.5);
            ylim(ylims);
            xlim(xlims1);
            box off;
            
            % Plot normalized bessel functions
            if is_normbessel
                figure(4); clf;
                b_dat = smooth(real(data1.coh_sum_win(ind)/data1.coh_num),npts_smooth);
                plot(faxis(ind),b_dat./abs(hilbert(b_dat)),'k','linewidth',3); hold on;
                plot(waxis/2/pi,b./SmoothAnalyticEnv(waxis/2/pi,b),'-r','linewidth',2); hold on;
                xlim(xlims);
                xlims1 = get(gca,'XLim');
                ylabel('J_{0}','fontsize',16);
                ylim([-2 2]);
            end
            
            if isfigure2
            f12 = figure(12);
            clf
            T = length(data1.coh_sum);
            faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(faxis>0);
            plot(faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),npts_smooth),'-r');
            xlim([frange_fit(1) frange_fit(2)])
            end
            psfile = [XSP_fig_path,'Xsp_',comp{1}(1),'_',sta1,'_',sta2,'_J0J1.pdf'];
            %print('-dpsc2',psfile);
            drawnow
            if isoutput
                save2pdf(psfile,f3,250);
            end

            
%             pause;
        end
        if isoutput
            save(sprintf('%s/%s_%s_xsp.mat',XSP_path,sta1,sta2),'xspinfo','twloc','waxis');
        end
        
%         pause;
        
        
    end %end of station j
end  %end of station i
