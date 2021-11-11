% Extract phase velocity dispersion between station pairs by fitting J0 bessel 
% function to real(ccf)
% Uses cross spectral fitting technique of Menke & Jin (2015) BSSA 
% DOI:10.1785/0120140245
%
% Define own starting phase velocity dispersion c manually or using 
% functions/calc_Rayleigh_disp for a simple layered model (does not work for 
% models with a water column).
%
% https://github.com/jbrussell
clear
close all;

setup_parameters;

%======================= PARAMETERS =======================%

% % RAYLEIGH FUND MODE (Mantle)
% comp = {'ZZ'}; %'RR'; 'ZZ'; 'TT'
% windir = 'window3hr';
% xspdir = 'ZZ_0S_LRT'; % output directory of phase velocities
% N_wl = 1; % Number of wavelengths required
% Ninterp = 25; % [] or Number of points to interpolate to;
% % Use picks from Linear Radon Transform? (./mat-LRTdisp/)
% is_LRT_picks = 1; % Use picks from Radon Transform to determine starting dispersion model and frequencies
% LRT_method = 'CGG_weight';
% mode_br = 0; % desired mode branch [0=fund.]
% xlims = [1/70 1/5]; % limits for plotting
% frange_LRT = [1/40 1/3]; % Frequency range of LRT panel for reading in picks
% frange_fit = [1/40 1/13]; % Frequency range to fit over! Can be more restrictive than where picks were made
% damp = [1; 1; 1]; % [fit, smoothness, slope]
% is_normbessel = 0; % normalize bessel function by analytic envelope?

% % RAYLEIGH FUND MODE (water and mantle)
% comp = {'PP'}; %'RR'; 'ZZ'; 'TT'
% windir = 'window3hr';
% xspdir = 'PP_0S_LRT'; % output directory of phase velocities
% N_wl = 1; % Number of wavelengths required
% Ninterp = 40; % [] or Number of points to interpolate to;
% % Use picks from Linear Radon Transform? (./mat-LRTdisp/)
% is_LRT_picks = 1; % Use picks from Radon Transform to determine starting dispersion model and frequencies
% LRT_method = 'CGG_weight';
% mode_br = 0; % desired mode branch [0=fund.]
% xlims = [1/70 1/5]; % limits for plotting
% frange_LRT = [1/40 1/3]; % Frequency range of LRT panel for reading in picks
% frange_fit = [1/40 1/5]; % Frequency range to fit over! Can be more restrictive than where picks were made
% damp = [1; 1; 1]; % [fit, smoothness, slope]
% is_normbessel = 0; % normalize bessel function by analytic envelope?

% RAYLEIGH 1st MODE (mantle)
comp = {'ZZ'}; %'RR'; 'ZZ'; 'TT'
windir = 'window3hr_Zcorr_tiltcomp';
xspdir = 'ZZ_1S_LRT'; % output directory of phase velocities
N_wl = 1; % Number of wavelengths required
Ninterp = 15; % [] or Number of points to interpolate to;
% Use picks from Linear Radon Transform? (./mat-LRTdisp/)
is_LRT_picks = 1; % Use picks from Radon Transform to determine starting dispersion model and frequencies
LRT_method = 'CGG_weight';
mode_br = 1; % desired mode branch [0=fund.]
xlims = [1/70 1/5]; % limits for plotting
frange_LRT = [1/40 1/3]; % Frequency range of LRT panel for reading in picks
frange_fit = [1/12 1/5]; % Frequency range to fit over! Can be more restrictive than where picks were made
damp = [1; 1; 1]; % [fit, smoothness, slope]
is_normbessel = 0; % normalize bessel function by analytic envelope?
snr_thresh = 10;

% % LOVE Fund MODE (mantle)
% comp = {'TT'}; %'RR'; 'ZZ'; 'TT'
% windir = 'window3hr';
% xspdir = 'TT_0T_LRT'; % output directory of phase velocities
% N_wl = 1; % Number of wavelengths required
% Ninterp = 15; % [] or Number of points to interpolate to;
% % Use picks from Linear Radon Transform? (./mat-LRTdisp/)
% is_LRT_picks = 1; % Use picks from Radon Transform to determine starting dispersion model and frequencies
% LRT_method = 'CGG_weight';
% mode_br = 0; % desired mode branch [0=fund.]
% xlims = [1/70 1/5]; % limits for plotting
% frange_LRT = [1/40 1/3]; % Frequency range of LRT panel for reading in picks
% frange_fit = [1/12 1/5]; % Frequency range to fit over! Can be more restrictive than where picks were made
% damp = [1; 1; 1]; % [fit, smoothness, slope]
% is_normbessel = 0; % normalize bessel function by analytic envelope?

if ~is_LRT_picks
    frange_fit = [1/40 1/10]; % frequency range over which to fit bessel function
%     xlims = [1/70 1/9];
    Npers = 18; % Number of periods
    t_vec_all = 1./flip(linspace(frange_fit(1) , frange_fit(2) ,Npers)); % periods at which to extract phase velocity
end

is_resume = 0; % Resume from last processed file (1) or overwrite (0)
iswin = 1; % Use the time-domain windowed ccfs?
npts_smooth = 1; % Smoothing of bessel function: 1 = no smoothing

isoutput = 1; % Save *.mat file with results?
nearstadist = 0;
IsFigure = 1;
isfigure2 = 0;
isfigure_snr = 0;

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

% figure output path
if iswin
    XSP_fig_path = ['./figs/',windir,'/fullStack/Xsp/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_',num2str(N_wl),'wl_',xspdir,'/BesselDist/'];
else
    XSP_fig_path = ['./figs/',windir,'/fullStack/Xsp/',num2str(1/frange_fit(2)),'_',num2str(1/frange_fit(1)),'s_',num2str(N_wl),'wl_',xspdir,'/BesselDist_nowin/'];
end

if ~exist(XSP_fig_path)
    mkdir(XSP_fig_path);
end




warning off; %#ok<WNOFF>

stalist = parameters.stalist;
nsta=parameters.nsta; % number of target stations to calculate for
Nfreq = length(t_vec_all);
mat = [];
for ip = 1:Nfreq
    mat(ip).period = t_vec_all(ip);
    mat(ip).omega = 2*pi/t_vec_all(ip);
    mat(ip).c_pre = c_all(ip);
    mat(ip).xsp_real = [];
    mat(ip).xsp_imag = [];
    mat(ip).r = [];
    mat(ip).azi = [];
end
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
        
%         % Check to see if we have already done this
%         if is_resume && exist([XSP_path,sta1,'_',sta2,'_xsp.mat'])
%             disp('Already fit this one!')
%             continue
%         end
        clear data1 xcorf1 xsp1 filename
        
        %%% --- Load in the ccf --- %%%
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        
        if ~exist(filename,'file')
            disp(['not exist ',filename])
            continue;
        end
        
        
        data1 = load(filename);
        r1 = distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2),referenceEllipsoid('GRS80'))/1000;
        az1 = azimuth(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2),referenceEllipsoid('GRS80'));
        groupv_max = data1.max_grv;
        groupv_min = data1.min_grv;
        
        if r1 < nearstadist
            continue;
        end
        
        % Index wavelength criterion
%         I_wl = r1 ./ (t_vec_all .* c_all) > N_wl;
        I_wl = true(size(t_vec_all));
        if sum(I_wl) <= 1
            I_wl(1) = 1;
            I_wl(2) = 1;
        end
        c = c_all(I_wl);
        t_vec = t_vec_all(I_wl);
        c_std = c_all_std(I_wl);
        
        wholesec = npts;
        
        
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
        xcorf1_imag = imag(xcorf(1:N));
        xcorf1_imag(1) = 0;
        
        if isfigure2 
            figure(1)
            T = length(xcorf1);
            dt = 1;
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
        faxis = [0:N-1]*1/wholesec;

%         tw1 = ones(1,tN)*r1./c;
        
        xsp_real_f = interp1(faxis*2*pi,xcorf1,2*pi./t_vec);
        xsp_real_f = smooth(xsp_real_f,npts_smooth);
        xsp_imag_f = interp1(faxis*2*pi,xcorf1_imag,2*pi./t_vec);
        xsp_imag_f = smooth(xsp_imag_f,npts_smooth);
        
        
        %% %%% Calculate SNR %%%
        xcorf1 = data1.coh_sum./data1.coh_num;
        xcorf1_filtered = tukey_filt( xcorf1,[min(t_vec) max(t_vec)],1,0.25 );
        [snr, signal_ind] = calc_SNR(xcorf1_filtered,groupv_min,groupv_max,r1,isfigure_snr);
        
        if snr < snr_thresh
            continue
        end
        %%
        
        for ip = 1:Nfreq
            mat(ip).xsp_real = [mat(ip).xsp_real; xsp_real_f(ip)];
            mat(ip).xsp_imag = [mat(ip).xsp_imag; xsp_imag_f(ip)];
            mat(ip).r = [mat(ip).r; r1];
            mat(ip).azi = [mat(ip).azi; az1];
        end
        
        %%
        
        if 0 % plot initial bessel
            tw_init = interp1(twloc,tw1(1:tN),waxis,'linear');
            x_init = waxis.*tw_init;
            A = 1;
            binit = besselj(0,x_init)*A;
            binit = binit./mean(abs(binit)).*mean([abs(xsp1)]);
            plot(waxis/2/pi,binit,'-k','linewidth',2); hold on;
        end

        
%         if isoutput
%             save2pdf(psfile,f3,250);
%         end

    end %end of station j
end  %end of station i

%% Fit bessel at each frequency of interest
figure(19); clf;
for ip = 1:length(mat)
    omega = mat(ip).omega;
    xsp_real = mat(ip).xsp_real;
    r = mat(ip).r;
    c_start = mat(ip).c_pre;
    for ii = 1:2
        tw = r./c_start;
        %%% - Invert for the bessel function 2x - %%%
        options = optimoptions(@lsqnonlin,'TolFun',1e-12,'MaxIter',1500,'MaxFunEvals',1500);
%         [tw,~,res,~,~,~,J] = lsqnonlin(@(x) besselerr_dist(x,[xsp_real],omega,r,is_normbessel),[tw],[tw]*0.8,[tw]*1.2,options);
%         c_start = r./tw;
        [c,~,res,~,~,~,J] = lsqnonlin(@(x) besselerr_dist(x,[xsp_real],omega,r,is_normbessel),[c_start],[c_start]*0.5,[c_start]*1.5,options);
        c_start = c;
    end
    % ESTIMATE ERROR BARS ON MODEL! :: JBR - 2/2020
    % Calculate data variance from residual following Menke QMDA book eq. (4.31)
    % s_d^2 = E / (N-M)
    sigma_d2 = res'*res / (length(res)-length(c));
    Cov_m = inv(J'*J)*sigma_d2;
    sigma_m_c = diag(Cov_m).^(1/2);
    
%     mat(ip).c = r./tw;
    mat(ip).c = c;
    mat(ip).c_std = sigma_m_c;
    
    plot(mat(ip).period,mat(ip).c_pre,'ok'); hold on;
    errorbar(mat(ip).period,mat(ip).c,mat(ip).c_std,'.r');
    xlabel('Period (s)');
    ylabel('Phase Velocity (km/s)');
end

%% Plot by frequency
figure(20); clf;

Nrow = 5;
Ncol = 5;
for ip = 1:length(mat)
    subplot(Nrow,Ncol,ip); hold on;
    
    c = mat(ip).c;
%     c_pre = 12;
    period = mat(ip).period;
    omega = mat(ip).omega;
    xsp_real = mat(ip).xsp_real;
    r = mat(ip).r;
    wavelength = c .* period;
    Nwl = r ./ wavelength;
    
    % Plot observed
    plot(r,xsp_real,'.');
    
    % Plot estimated
    bessel_est = besselj(0,omega.*r./c);
    A_fac = mean(abs(xsp_real)) ./mean(abs(bessel_est));
    bessel_est = bessel_est .* A_fac;
    plot(r,bessel_est,'.r');
    
    % Plot model
    rplot = [min(r):1:max(r)];
    Nwlplot = rplot ./ c ./ period;
    bessel_plot = A_fac .* besselj(0,omega.*rplot./c);
    plot(rplot,bessel_plot,'-k');
   
    
    title([num2str(period),' s'])
    xlabel('Distance (km)');
    ylabel('Re(\rho)');
end

%% Plot all on one plot
figure(21); clf;
hold on;
for ip = 1:length(mat)
    
    c = mat(ip).c;
%     c_pre = 12;
    period = mat(ip).period;
    omega = mat(ip).omega;
    xsp_real = mat(ip).xsp_real;
    r = mat(ip).r;
    wavelength = c .* period;
    Nwl = r ./ wavelength;
    
    % Plot observed
    plot(Nwl,xsp_real,'.b');
    
    % Plot estimated
    bessel_est = besselj(0,omega.*r./c);
    A_fac = mean(abs(xsp_real)) ./mean(abs(bessel_est));
    bessel_est = bessel_est .* A_fac;
    plot(Nwl,bessel_est,'.r');
    
    % Plot model
    rplot = [min(r):1:max(r)];
    Nwlplot = rplot ./ c ./ period;
    bessel_plot = A_fac .* besselj(0,omega.*rplot./c);
    plot(Nwlplot,bessel_plot,'-k');
   
    
    title([num2str(period),' s'])
    xlabel('Distance/Wavelength');
    ylabel('Re(\rho)');
end

%% Polar plots by frequency (REAL)
figure(22); clf;
set(gcf,'position',[59           4         837        1021]);

Nrow = 5;
Ncol = 5;
sgtitle('Real','fontsize',25,'fontweight','bold');
for ip = per_ind %1:2:length(mat)
    ax = subplot(Nrow,Ncol,ip); %hold on;
    
    c = mat(ip).c;
%     c_pre = 12;
    period = mat(ip).period;
    omega = mat(ip).omega;
    xsp_real = mat(ip).xsp_real;
    r = mat(ip).r;
    azi = mat(ip).azi;
    wavelength = c .* period;
    Nwl = r ./ wavelength;
    
    % Plot observed
%     plot(r,xsp_real,'.');
    polarscatter(azi,r,25,xsp_real,'filled'); hold on;
    polarscatter(azi+180,r,25,xsp_real,'filled'); hold on;
    colormap(jet); 
    cb = colorbar;
    ylabel(cb,'Real Cross-spectrum');
    caxis(max(max([mat(:).xsp_real]))*0.5*[-1 1]);
    rlim([0 max(max([mat(:).r]))]);
    title([num2str(round(mat(ip).period)),' s'])
    set(gca,'ThetaZeroLocation','top','RTickLabel',[],'fontsize',15,'linewidth',2);
end

save2pdf(['./figs_paper/polarNCF_imag.pdf'],22,250);

%% Polar plots by frequency (IMAGINARY)
figure(23); clf;
set(gcf,'position',[59           4         837        1021]);

Nrow = 5;
Ncol = 5;
sgtitle('Imaginary','fontsize',25,'fontweight','bold');
for ip = per_ind %1:2:length(mat)
    ax = subplot(Nrow,Ncol,ip); %hold on;
    
    c = mat(ip).c;
%     c_pre = 12;
    period = mat(ip).period;
    omega = mat(ip).omega;
    xsp_imag = mat(ip).xsp_imag;
    r = mat(ip).r;
    azi = mat(ip).azi;
    wavelength = c .* period;
    Nwl = r ./ wavelength;
    
    % Plot observed
%     plot(r,xsp_real,'.');
    polarscatter(azi,r,25,xsp_imag,'filled'); hold on;
    polarscatter(azi+180,r,25,xsp_imag,'filled'); hold on;
    colormap(jet); 
    cb = colorbar;
    ylabel(cb,'Imaginary Cross-spectrum');
    caxis(max(max([mat(:).xsp_real]))*0.5*[-1 1]);
    rlim([0 max(max([mat(:).r]))]);
    title([num2str(round(mat(ip).period)),' s'])
    set(gca,'ThetaZeroLocation','top','RTickLabel',[],'fontsize',15,'linewidth',2);
end

save2pdf(['./figs_paper/polarNCF_imag.pdf'],23,250);
