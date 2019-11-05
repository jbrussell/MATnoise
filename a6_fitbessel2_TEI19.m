% Program to fit the xcorf and get the travel time for each frequency for two station pairs
% 
% JBR 5/19/18 : This version fits J0 to the real part and J-1 to the
% imaginary part. The starting dispersion model is taken from a given
%

clear
close all;

global tN
global waxis
global twloc
global weight


setup_parameters;

% % FUNDAMENTAL MODE RAYLEIGH (NARROW WINDOW)
% comp = {'ZZ'}; %'RR'; 'ZZ'; 'TT'
% windir = 'window3hr_Z'; %'window3hr_LH_Zcorr'; %'window3hr_LH_Zcorr'; %'window0.2hr'; %'window24hr_specwhite';
% xspdir = 'test_1.6win'; %'test_1.6win'; %'test_1.6nowin'; %'test_1.6win';
% % qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3_h2o4.5km.s0to200.q';
% qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3_h2o4.5km.s0to250.q';
% mode = 0; % 0 = fundamental; % of starting phv
% frange = [1/25 1/14]; %[1/35 1/14]; %[0.1 0.25];
% % For SNR calculation
% groupv_max = inf;
% groupv_min = 1.6;
% Npers = 12;
% xlims = [0.02 0.1];

% % FUNDAMENTAL MODE RAYLEIGH (WIDE WINDOW)
% comp = {'ZZ'}; %'RR'; 'ZZ'; 'TT'
% windir = 'window3hr_0.8kms'; %'window3hr_LH_Zcorr'; %'window3hr_LH_Zcorr'; %'window0.2hr'; %'window24hr_specwhite';
% xspdir = 'test_0.8win'; %'test_1.6win'; %'test_1.6nowin'; %'test_1.6win';
% % qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3.s0to200.q'; % for starting phv
% qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3_h2o4.5km.s0to200.q';
% mode = 0; % 0 = fundamental; % of starting phv
% frange = [1/25 1/14]; %[1/35 1/14]; %[0.1 0.25];
% % For SNR calculation
% groupv_max = inf;
% groupv_min = 0.8;
% Npers = 12;
% xlims = [0.02 0.1];

% % FIRST OVERTONE RAYLEIGH
% comp = {'ZZ'}; %'RR'; 'ZZ'; 'TT'
% windir = 'window3hr_Zcorr_tiltcomp'; %'window3hr_LH_Zcorr'; %'window3hr_LH_Zcorr'; %'window0.2hr'; %'window24hr_specwhite';
% xspdir = 'test_1.4win'; %'test_1.6win'; %'test_1.6nowin'; %'test_1.6win';
% % qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3.s0to200.q'; % for starting phv
% % qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3_h2o4.5km.s0to200.q';
% qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3_h2o4.5km.s0to250.q';
% mode = 1; % 0 = fundamental; % of starting phv
% frange = [1/10 1/5]; %[0.1 0.25];
% % For SNR calculation
% groupv_max = inf;
% groupv_min = 1.4;
% Npers = 12;
% xlims = [1/12 1/3];

% % FUNDAMENTAL MODE LOVE
% comp = {'TT'}; %'RR'; 'ZZ'; 'TT'
% windir = 'window3hr'; %'window3hr_LH_Zcorr'; %'window3hr_LH_Zcorr'; %'window0.2hr'; %'window24hr_specwhite';
% xspdir = 'test_1.4win'; %'test_1.6win'; %'test_1.6nowin'; %'test_1.6win';
% % qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3.s0to200.q'; % for starting phv
% % qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3_h2o4.5km.s0to200.q';
% qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3.t0to333.q';
% mode = 0; % 0 = fundamental; % of starting phv
% frange = [1/9 1/3]; %[0.1 0.25];
% % For SNR calculation
% groupv_max = inf;
% groupv_min = 1.4;
% Npers = 21;
% xlims = [1/12 1/3];
% t_vec = 1./flip(linspace(frange(1) , frange(2) ,Npers));
% % c = [3.7747 % fast starting c
% %     3.8080
% %     3.8421
% %     3.8772
% %     3.9126
% %     3.9484
% %     3.9846
% %     4.0217
% %     4.0599
% %     4.0991
% %     4.1394
% %     4.1807
% %     4.2232
% %     4.2668
% %     4.3117
% %     4.3577
% %     4.4051
% %     4.4537
% %     4.5032
% %     4.5534
% %     4.6045]';
% c = [3.6228 % slow starting c
%     3.6582
%     3.6948
%     3.7327
%     3.7709
%     3.8094
%     3.8482
%     3.8883
%     3.9298
%     3.9727
%     4.0170
%     4.0629
%     4.1103
%     4.1594
%     4.2102
%     4.2627
%     4.3171
%     4.3734
%     4.4307
%     4.4889
%     4.5482]';

% FUNDAMENTAL MODE LOVE (mean starting c)
comp = {'TT'}; %'RR'; 'ZZ'; 'TT'
windir = 'window3hr'; %'window3hr_LH_Zcorr'; %'window3hr_LH_Zcorr'; %'window0.2hr'; %'window24hr_specwhite';
xspdir = 'test_1.4win_meanc'; %'test_1.6win'; %'test_1.6nowin'; %'test_1.6win';
% qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3.s0to200.q'; % for starting phv
% qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3_h2o4.5km.s0to200.q';
qfile = 'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3.t0to333.q';
mode = 0; % 0 = fundamental; % of starting phv
frange = [1/9 1/3]; %[0.1 0.25];
% For SNR calculation
groupv_max = inf;
groupv_min = 1.4;
Npers = 21;
xlims = [1/12 1/3];
t_vec = 1./flip(linspace(frange(1) , frange(2) ,Npers));
c = [3.6639
    3.6953
    3.7284
    3.7632
    3.7976
    3.8316
    3.8654
    3.9010
    3.9387
    3.9784
    4.0200
    4.0636
    4.1092
    4.1569
    4.2067
    4.2588
    4.3131
    4.3698
    4.4269
    4.4844
    4.5422]';

iswin = 1; % Use the time-domain windowed ccfs?
npts_smooth = 1; % 1 = no smoothing

isoutput = 1; % Save *.mat file with results?
nearstadist = 0;
IsFigure = 1;
isplotinit = 0;
isfigure2 = 0;
isfigure_snr = 1;

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

%% Make an initial model
% From GOC_CC card
vec_h = [3 2 4 12]; % Layer thickness
vec_vs = [1.1 1.2 2.8 3.7 4.6];
vec_vp = vec_vs.*1.8; vec_pv(1) = 1.5;
vec_rho = [1.03 1.5 3.02 3.027 3.342];
vr = mat_disperse(vec_h,vec_rho,vec_vp,vec_vs,1./t_vec);
c = vr(:,1)';
c_start = c;

% c = [3.4837    3.6341    3.7458    3.8223    3.8878    3.9451   4.0013    4.0522    4.0951    4.1337    4.1683    4.2098]; % 'test_1.6win' avg
% 10.0000   11.2063   12.5580   14.0729   15.7704   17.6727   19.8045   22.1935   24.8706   27.8706   31.2325   35.0000
% c = [3    3.2    3.3    3.4464    3.8400    3.9589    4.0097    4.0363    4.0515    4.0600    4.0644    4.0661];

if exist('c','var') == 0 % check if phase velocities exist, if not read them in
    [~,~,c] = readMINEOS_qfile2(qfile,t_vec,mode);
end
c_start = c;

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
%         snrdata = real(ifft((xcorf1_filtered)));
%         snrdata = fftshift(snrdata);
%         groupv_max = 10;
%         groupv_min = 1.6;
%         NN= length(snrdata);
%         lag = [-floor(NN/2):floor(NN/2)];
%         win_min = r1./groupv_max;
%         win_max = r1./groupv_min;
%         if win_min < 15; win_min = 0; end
%         if win_max < 50; win_max = 50; end
%         signal_ind = ((lag>-win_max & lag<-win_min) | (lag>win_min & lag<win_max));
%         signal_amp = sum(snrdata(signal_ind).^2)/length(snrdata(signal_ind));
%         noise_amp = sum(snrdata(~signal_ind).^2)/length(snrdata(~signal_ind));
%         snr = signal_amp/noise_amp;
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
                %plot(wvec1/2/pi,binit,'-k')
                plot(waxis/2/pi,binit,'-k','linewidth',2); hold on;
            end
            
            f3 = figure(3); clf; hold on; 
%             set(gcf, 'Color', 'w','position',[259    63   713   642]);
%             set(gcf,'color','w','Position',[289   218   532   487]);
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
%             xlim([frange(1) frange(2)])
            xlim(xlims);
            xlims1 = get(gca,'XLim');
%             xlabel('Frequency (Hz)','fontsize',16);
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
%stapairn
%%soundsc(rand(2000,1),1000,8)
