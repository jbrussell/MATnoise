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

global tN
global waxis
global twloc
global weight
setup_parameters;

%======================= PARAMETERS =======================%
comp = {'ZZ'}; %'RR'; 'ZZ'; 'TT'
windir = 'window3hr';
xspdir = 'phv_dir'; % output directory of phase velocities
frange = [1/5 1/10]; % frequency range over which to fit bessel function

% Windowing uptions for SNR calculation
groupv_max = inf;
groupv_min = 1.4;

Npers = 21; % Number of periods
xlims = [1/12 1/3];
t_vec = 1./flip(linspace(frange(1) , frange(2) ,Npers)); % periods at which to extract phase velocity

% % Manually input phase velocity
% c = [3.6639
%     3.6953
%     3.7284
%     3.7632
%     3.7976
%     3.8316
%     3.8654
%     3.9010
%     3.9387
%     3.9784
%     4.0200
%     4.0636
%     4.1092
%     4.1569
%     4.2067
%     4.2588
%     4.3131
%     4.3698
%     4.4269
%     4.4844
%     4.5422]';

iswin = 1; % Use the time-domain windowed ccfs?
npts_smooth = 1; % 1 = no smoothing

isoutput = 1; % Save *.mat file with results?
nearstadist = 0;
IsFigure = 1;
isplotinit = 0;
isfigure2 = 0;
isfigure_snr = 1;

%% Make the initial phase velocity dispersion model

% calc_Rayleigh_disp
vec_h = [3 2 4 12]; % Layer thickness
vec_vs = [1.1 1.2 2.8 3.7 4.6];
vec_vp = vec_vs.*1.8; vec_pv(1) = 1.5;
vec_rho = [1.03 1.5 3.02 3.027 3.342];
vr = mat_disperse(vec_h,vec_rho,vec_vp,vec_vs,1./t_vec);
c = vr(:,1)';
c_start = c;

% Manually
% c = [3.4837    3.6341    3.7458    3.8223    3.8878    3.9451   4.0013    4.0522    4.0951    4.1337    4.1683    4.2098]; % 'test_1.6win' avg
% 10.0000   11.2063   12.5580   14.0729   15.7704   17.6727   19.8045   22.1935   24.8706   27.8706   31.2325   35.0000
% c = [3    3.2    3.3    3.4464    3.8400    3.9589    4.0097    4.0363    4.0515    4.0600    4.0644    4.0661];

% From MINEOS .q file (https://github.com/jbrussell/MINEOS_synthetics)
% if exist('c','var') == 0 % check if phase velocities exist, if not read them in
%     [~,~,c] = readMINEOS_qfile2(qfile,t_vec,mode);
% end
% c_start = c;
%==========================================================%

twloc=1./t_vec;

% LIMITS
if comp{1}(1) == 'R'
    ylims = [3.2 4.5];
elseif comp{1}(1) == 'Z'
%     ylims = [1.5 5.0];
    ylims = [2 4.5];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tN = length(t_vec);
wholesec = npts;
wvec1 = (2*pi)./t_vec;
wvec1 = wvec1';
refc = c';

% input path
%ccf_path = ['./ccf/',windir,'/fullStack/ccf',comp{1},'/'];
ccf_path = [parameters.ccfpath,windir,'/fullStack/ccf',comp{1},'/'];

% output path
XSP_path = ['./Xsp/',windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/'];

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
    XSP_fig_path = ['./figs/',windir,'/fullStack/Xsp/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/TEI19/'];
else
    XSP_fig_path = ['./figs/',windir,'/fullStack/Xsp/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/TEI19_nowin/'];
end

if ~exist(XSP_fig_path)
    mkdir(XSP_fig_path);
end




warning off; %#ok<WNOFF>




% Get your axis correct
twloc = twloc*2*pi;
waxis = (frange(1):1/wholesec:frange(2))*2*pi;

%faxis = [0:wholesec/2 -1]*1/wholesec;

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
        if 0%exist([XSP_path,sta1,'_',sta2,'_xsp.mat'])
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
        delta = distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2));
        r1    = deg2km(delta); % distance
        
        if r1 < nearstadist
            continue;
        end
        
        
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
            dt = 1;
            temp_faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(temp_faxis>0);
            subplot(2,1,1)
            %plot(temp_faxis(ind),smooth(real(xcorf1(ind)),100));
            plot(flip(temp_faxis(ind),smooth(real(xcorf1(ind)),50)));
            xlim([frange(1) frange(2)])
            hold on
            subplot(2,1,2)
            %plot(temp_faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),100),'-r')
            plot(temp_faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),50),'-r');
            xlim([frange(1) frange(2)])
            
        end

        %%% - Convert xcorf into spherical frequency - %%%
        faxis = [0:N-1]*1/wholesec;
        xsp1 = interp1(faxis*2*pi,xcorf1,waxis);

        %xsp1 = smooth(xsp1,50);
        xsp1 = smooth(xsp1,npts_smooth);

        tw1 = ones(1,tN)*r1./c;
        
        %%% - Invert for the bessel function 2x - %%%
        options = optimoptions(@lsqnonlin,'TolFun',1e-12,'MaxIter',1500,'MaxFunEvals',1500);
        weight  = 1./waxis;
        tw2 = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw1],[tw1]*0.8,[tw1]*1.2,options);
%         tw2 = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw1],[],[],options);
        
        weight(:) = 1;
        tw = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw2],[tw2]*0.8,[tw2]*1.2,options);
%         tw = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw2],[],[],options);
        
        
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
        xspinfo.coherenum = data1.coh_num;
        err = besselerr(tw,xsp1);
        err = err(1:length(waxis));
        xspinfo.sumerr = sum(err.^2)./sum((xsp1./weight(:)).^2);
        xspinfo.err = err./weight(:);
        xspinfo.tw1 = tw1;
        xspinfo.twloc = twloc;
        xspinfo.c = r1./tw;
        xspinfo.per = 1./(twloc/2/pi);
        xspinfo.c_start = c_start;
        
        data = r1./tw;
        

        %% %%% Calculate SNR %%%
        xcorf1 = data1.coh_sum./data1.coh_num;
        xcorf1_filtered = tukey_filt( xcorf1,flip(1./frange),1,0.25 );
        [snr, signal_ind] = calc_SNR(xcorf1_filtered,groupv_min,groupv_max,r1,isfigure_snr);
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
            plot_SNR(xcorf1_filtered,groupv_min,groupv_max,r1,ax1);
            set(ax1,'box','off');

            % REAL PART (J0)
            subplot(3,1,2); box off;
            tww = interp1(twloc,tw(1:tN),waxis,'linear');
            x = waxis.*tww;
            A = 1;
            b = besselj(0,x)*A;
            b = b./mean(abs(b)).*mean([abs(xsp1)]);           
            T = length(data1.coh_sum);
            dt = 1;
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
            plot(twloc/2/pi,r1./tw1,'o-','color',[0.5 0.5 0.5],'linewidth',2);hold on;
            plot(twloc/2/pi,r1./tw,'o-','color',[1 0 0 ],'linewidth',2);
%             errorbar(twloc/2/pi,r1./tw,xspinfo.err,'ro-','linewidth',2);
            title([sta1,'-',sta2],'fontsize',16)
            xlabel('Frequency (Hz)','fontsize',16);
            ylabel('Phase Velocity (km/s)','fontsize',16);
            set(gca,'fontsize',16,'linewidth',1.5);
            ylim(ylims);
            xlim(xlims1);
            box off;
            
            if isfigure2
            f12 = figure(12);
            clf
            T = length(data1.coh_sum);
            dt = 1;
            faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(faxis>0);
            plot(faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),npts_smooth),'-r');
            xlim([frange(1) frange(2)])
            end
            psfile = [XSP_fig_path,'Xsp_',comp{1}(1),'_',sta1,'_',sta2,'_J0J1.pdf'];
            %print('-dpsc2',psfile);
            save2pdf(psfile,f3,1000);

            
%             pause;
        end
        if isoutput
            save(sprintf('%s/%s_%s_xsp.mat',XSP_path,sta1,sta2),'xspinfo','twloc','waxis');
        end
        
%         pause;
        
        
    end %end of station j
end  %end of station i
