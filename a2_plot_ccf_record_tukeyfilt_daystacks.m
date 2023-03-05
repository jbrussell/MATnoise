% Plot the daily cross-spectra in the time domain for individual station
% pairs
%
% JBR, 3/5/2023
clear all;
setup_parameters;
IsFigure = 0;
IsFigure_GAUS = 0; % Plot frequency domain filtered and unfiltered

%======================= PARAMETERS =======================%
comp = 'ZZ'; %'ZZ'; %'RR'; %'TT';
coperiod = [1 10]; %[4 10]; %[5 25]; %[2 8]; %[2 16]; %[5 25]; %[ 8 30 ]; % Periods to filter between

%%% --- Parameters to build up gaussian filters --- %%% 
% (effects the width of the filter in the frequency domain)
costap_wid = 0.25; % 0 => box filter; 1 => Hann window

% Window Velocities
max_grv = 5.5; %8.0;
min_grv = 2.2; %1.6; % FOR WINDOWING!

%==========================================================%
xlims = [-500 500];

dt = parameters.dt;
stalist = parameters.stalist;
% stalist =  stalist(1:21);
nsta = length(stalist);
winlength = parameters.winlength;
figpath = parameters.figpath;
%fig_winlength_path = [figpath,'window',num2str(winlength),'hr/fullStack/'];
% custom directory names
    windir = 'window3hr'; %'window0.2hr'; %'window24hr_specwhite';
    fig_winlength_path = [figpath,windir,'/dayStack/'];

%------------ PATH INFORMATION -------------%
% ccf_path = parameters.ccfpath;
ccf_path = './ccf/';
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

ccf_path = [ccf_stack_path,'ccf',comp,'/',];
npairall = 0;
sta1sta2_dist = [];
%------------ LOAD DATA AND PLOT IN TIME DOMAIN -------------%
for ista1=1:nsta % loop over all stations
    sta1=char(stalist(ista1,:));
    
    for ista2 = 1: nsta % loop over station pairs
        sta2 = char(stalist(ista2,:));
        
        npairall = 0;
        ccf_stack = {};
        nstapair = 0;
    %         sta2 = char(stalist(ista2,:));

        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        
        ccf_day_path = [ccf_daystack_path,'ccf',comp,'/',sta1,'/',sta1,'_',sta2,'_*.mat'];
        fils = dir(ccf_day_path);
        if isempty(fils)
            disp(['No files for ',sta1,'-',sta2,'... skipping'])
            continue
        end
        
        nstapair = nstapair + 1;
        npairall = npairall + 1; % number of total station pairs

        % Get correct station lats & lons
        Ista1 = find(strcmp(stalist,sta1));
        Ista2 = find(strcmp(stalist,sta2));
        sta1sta2_dist(nstapair) = deg2km(distance(stalat(Ista1),stalon(Ista1),stalat(Ista2),stalon(Ista2)));
        
        
        din = [];
        timeflag = [];
        tvec = NaT(length(fils),1);
        daySNR = [];
    %     sta1sta2_dist = [];
        for ifil = 1:length(fils)
            fldr = [ccf_daystack_path,'ccf',comp,'/',sta1,'/'];
            temp = load([fldr,'/',fils(ifil).name]);
        %     temp = load([fils(ifil).folder,'/',fils(ifil).name]);
            tvec(ifil) = temp.starttime;
            ccf_day = temp.coh_sum./temp.coh_num_day;
            [ ccf_filtered ] = tukey_filt( ccf_day,coperiod,dt,costap_wid );

    %             % Do butterworth filtering after rearranging
    %             [b, a] = butter(4,flip(1./coperiod)*2*dt); % Butterworth Filter
    %             ccf_ifft = filtfilt(b,a,real(ifft(ccf_day)));
    %             ccf_filtered = fft(ccf_ifft);

            ccf_day = ccf_filtered;
            N = length(ccf_day);
            ccf_day_ifft = real(ifft(2*ccf_day([1:N/2+1]),N)); % inverse FFT to get time domain
            %rearrange and keep values corresponding to lags: -(len-1):+(len-1)
            ccf_day_ifft = [ccf_day_ifft(end-N+2:end)' ; ccf_day_ifft(1:N)'];

            % Build time axis
            N = length(ccf_day_ifft);
            time = ([0:N-1]-floor(N/2))*dt;
            timeall = [time(time<0), time(time>=0)];
            Ikeep = abs(timeall)<=500; % index data to keep
            timeflag(:,ifil) = timeall(Ikeep);
            din(:,ifil) = ccf_day_ifft(Ikeep);

            % Calculate SNR
            r = sta1sta2_dist(nstapair);
            tmin = -r/min_grv; %50;
            tmax = r/min_grv; %150;
            Isignal = timeflag(:,ifil)>=tmin & timeflag(:,ifil)<=tmax;
            signal = din(Isignal,ifil);
            noise = din(~Isignal,ifil);
            daySNR(ifil) = max(abs(signal))/sqrt(mean(noise.^2));
    %             daySNR(ifil) = sqrt(mean(signal.^2))/sqrt(mean(noise.^2));
        end
        dayvec = [0; cumsum(seconds(diff(tvec))/60/60/24)]'; %[1:Ndays];
        % Calculate reference stack
        dstack = sum(din(:,:),2)./max(abs(sum(din(:,:),2)));
        ccf_stack{npairall} = dstack;

        twin = r/min_grv;
        if twin < 50
            twin = 50;
        end
        Igrvwin = timeflag(:,1)<=twin & timeflag(:,1)>=-twin;
        Igrvwinmat = repmat(Igrvwin,1,size(din,2));
        norm_mat = repmat(max(abs(din.*Igrvwinmat)),size(din,1),1);
        din_norm = din./norm_mat; din_norm(isnan(din_norm))=0;
    %         [drobust_norm,~]=smartstack(din_norm,'stackwindow',[robustitmin,robustitmax],'stacktype','robust');
    %         ccf_smartstack_all_norm{npairall} = drobust_norm/max(abs(drobust_norm));

        f101 = figure(101); clf; hold on;
        set(gcf,'color','w','position',[340   164   572   995]);
        N= length(ccf_stack{npairall});
        % time = [-N/2:N/2]*dt; % build lagtime vector for plotting
        time = ([0:N-1]-floor(N/2))*dt;
        indtime = find(abs(time)<=500);
        imagesc(timeflag(:,1)',dayvec,din_norm');
    %         imagesc(timeflag(:,1)',[1:length(fils)],din');
    %         plot(time(indtime),ccf_smartstack_all{npairall}(indtime)*10 - 10,'-r','linewidth',1);
        plot(time(indtime),dstack(indtime)*10-10,'-k','linewidth',1);
        xlim([-250 250]);
        ylim([-20 max(dayvec)]);
        xlabel('Lag time (s)','fontsize',15);
        ylabel('Day of deployment','fontsize',15);
        title([sta1,'-',sta2,' ',num2str(round(r)),' km (',num2str(coperiod(1)),'-',num2str(coperiod(2)),'s)'],'fontweight','bold','fontsize',15);
        colormap(redbluecmap);
        set(gca,'fontsize',15,'linewidth',1.5,'tickDir','out','box','on','layer','top','ydir','reverse')
        caxis([-1 1]);
    %         datetick('y','keepticks','keeplimits')
        save2pdf([figpath,'ccf',comp,'_',sta1,'_',sta2,'.pdf'],f101,100);
    end % ista2
end % ista1



