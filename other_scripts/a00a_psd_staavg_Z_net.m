% Script to calculate average station spectra so we can compare values across the
% whole array
%
% jbrussell - 11/2025

clear;
setup_parameters;

strSACcomp = 'BHZ';
strNAMEcomp = 'ZZ';
IsFigure1 = 1;
IsFigure2 = 0;

% OUTPUT SETTINGS
IsOutputFullstack = 1; % Save full year ccf stacks
IsOutputMonthstack = 1; % save month ccf stacks
IsOutputDaystack = 0; % save day ccf stacks
IsOutputSinglestack = 0; % save single ccf before stacking
IsOutputSeismograms = 0; % save raw seismograms before cross-correlating

% GENERAL PROCESSING
IsRemoveIR = 0; % remove instrument response
units_RemoveIR = 'M'; % 'M' displacement | 'M/S' velocity
IsDetrend = 1; % detrend the data
IsTaper = 1; % Apply cosine taper to data chunks

%%%%%%%%%%% OPTIONS FOR PREPROCESSING %%%%%%%%%%%%
% (1) ONE-BIT NORMALIZATION & SPECTRAL WHITENING? (Bensen et al. 2007)
IsSpecWhiten = 0; % Whiten spectrum
IsOBN = 0; % One-bit normalization

% (2) TIME-FREQUENCY NORMALIZATION (Ekstrom et al. 2009; Shen et al. 2011)
IsFTN = 0; % Frequency-time normalization? (If 1, applied instead of whitening and one-bit normalization)
frange_FTN = [1/60 1/10]; % frequency range over which to construct FTN seismograms

% (3) BASIC PREFILTER (Ekstrom 2011)
IsPrefilter = 0; % apply butterworth bandpass filter before cross-correlation?
frange_prefilt = [1/100 1/10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Setup parallel pool
% Nworkers = 4; % number of workers in pool for parallel processing
delete(gcp('nocreate'));
% parpool(Nworkers);

% input path
datadir = parameters.datapath;
PZpath = parameters.PZpath;
figpath = [parameters.figpath,'/psd/'];
seis_path = parameters.seis_path;
orientation_path = parameters.orientation_path;
dt = parameters.dt;
winlength = parameters.winlength;

year = ''; %'2012';
Nstart_sec = parameters.Nstart_sec; % (seconds) offset start of file
Nstart = Nstart_sec/dt; % Number of samples
comp = parameters.comp;

%dist_min = 20;
dist_min = parameters.mindist;

% Build File Structure: cross-correlations
% psd_path = parameters.ccfpath;
psd_path = ['./psd/'];
psd_winlength_path = [psd_path,'window',num2str(winlength),'hr/'];
psd_singlestack_path = [psd_winlength_path,'single/'];
psd_singlestack_path = [psd_winlength_path,'dayStack/'];
psd_monthstack_path = [psd_winlength_path,'monthStack/'];
psd_fullstack_path = [psd_winlength_path,'fullStack/'];

if ~exist(psd_path)
    mkdir(psd_path)
end
if ~exist(psd_winlength_path)
    mkdir(psd_winlength_path)
end
if ~exist(psd_singlestack_path)
    mkdir(psd_singlestack_path)
end
if ~exist(psd_singlestack_path)
    mkdir(psd_singlestack_path)
end
if ~exist(psd_monthstack_path)
    mkdir(psd_monthstack_path)
end
if ~exist(psd_fullstack_path)
    mkdir(psd_fullstack_path)
end

PATHS = {psd_singlestack_path; psd_singlestack_path; psd_monthstack_path; psd_fullstack_path};
for ipath = 1:length(PATHS)
    psdZ_path = [PATHS{ipath},'psd',strNAMEcomp,'/'];
    if ~exist(psdZ_path)
        mkdir(psdZ_path);
    end
end

% Build File Structure: figures
fig_winlength_path = [figpath,'window',num2str(winlength),'hr/'];
if ~exist(figpath)
    mkdir(figpath);
end
if ~exist(fig_winlength_path)
    mkdir(fig_winlength_path)
end

% Build File Structure: windowed seismograms
seis_winlength_path = [seis_path,'window',num2str(winlength),'hr/'];
if ~exist(seis_path)
    mkdir(seis_path);
end
if ~exist(seis_winlength_path)
    mkdir(seis_winlength_path)
end

warning('off','MATLAB:nargchk:deprecated')
%% ------------------- loop through center station station-------------------

stalist = parameters.stalist;
netlist = parameters.netlist;
nsta=parameters.nsta; % number of target stations to calculate for

% READ OBS ORIENTATIONS
% [slist, orientations] = textread(orientation_path,'%s%f\n');

% Calculate filter coefficients for FTN
if IsFTN
    [ b, a ] = get_filter_TFcoeffs( frange_FTN, dt );
end

% Get list of stations that have already been calculated (the last one in
% the list will still be in progress)

files = dir([psd_fullstack_path,'psd',strNAMEcomp,'/']);
for ii=3:length(files)
    check = dir([psd_fullstack_path,'psd',strNAMEcomp,'/',files(ii).name,'/*.mat']);
    if length(check)==0
        rmdir([psd_fullstack_path,'psd',strNAMEcomp,'/',files(ii).name])
        display(['removing empty folders ',files(ii).name])
    end
end
files = dir([psd_fullstack_path,'psd',strNAMEcomp,'/']);
sta_processed = {files([files.isdir]).name}; % gather folder names
sta_processed = sta_processed(~ismember(sta_processed ,{'.','..'})); % remove annoying '.' and '..' entries
sta_processed = sta_processed(1:end-1); % get rid of last one from list since it's still in progress

% Delete TA.H17A and (FX...) from sta_processed and rerun

for ista1=1:nsta

    sta1=char(stalist(ista1,:));
    net1=char(netlist(ista1,:));
    toks = strsplit(sta1,'.');
    stacode1 = toks{2};

    % Check if station has already been processed. If so, skip
    if any(strcmp(sta1,sta_processed))
        disp([sta1,' already processed. Skipping...'])
        continue
    end

    % Build station directories
    for ipath = 1:length(PATHS)
        psdZ_path = [PATHS{ipath},'psd',strNAMEcomp,'/'];
        if ~exist([psdZ_path,sta1])
            mkdir([psdZ_path,sta1]);
        end
    end
    seisZ_path = [seis_winlength_path,strNAMEcomp(1),'/'];
    if ~exist([seisZ_path,sta1])
        mkdir([seisZ_path,sta1]);
    end

    list1 = dir([datadir,'/',net1,'/',stacode1,'/*',strSACcomp,'.sac']);

    if isempty(list1)
        display(['No data for ',sta1]);
        continue
    end

    clear lat1 lat2 lon1 lon2 dist az baz vec_tz2 Z2raw vec_tz Z1raw



    % check to see if we've already done this ccf
    if exist([psdZ_path,sta1,'/',sta1,'_f.mat'])
        display(['PSD already exists for ',sta1,'... skip this station']);
        continue
    end

    display(['Calculating PSDs for staion : ',sta1]);
    % -------------loop through each half day--------------------
    nday_stack=0;
    fftSZ = 0;
    coh_num = 0;

    % Get a list of all available data
    ihday = 0;
    ii_fftday = 0;
    month_counter = 0;
    imonth = 0;
    ii_fftday_month = 0;
    clear fftS1Z_day_sum
    for ifil = 1:length(list1)
        file1cZ = list1(ifil).name;

        str = strsplit(file1cZ,'.');
        hdayid = [str{2},'.',str{3},'.',str{4},'.',str{5},'.',str{6}];
      
        if month_counter == 0
            fftSZ_month = 0;
            coh_num_month = 0;
        end
        clear data1cZ 
        ihday = ihday +1;
        month_counter = month_counter + 1;
        clear temp
        %temp = strsplit(daylist1(ihday).name,'.');

        disp(['Looking at ',hdayid,' ',sta1]);

        data1cZ= dir([datadir,'/',net1,'/',stacode1,'/',year,'/',stacode1,'.',hdayid,'.*',strSACcomp,'.sac']);

        data1cZ =  [datadir,'/',net1,'/',stacode1,'/',year,'/',data1cZ.name];

        %------------------- TEST IF DATA EXIST------------------------
        try
            [S1Zt,S1Zraw,S1,S1Ztstart] = load_sac(data1cZ);
        catch
            disp(['Something wrong sac file... skipping ',data1cZ])
            continue
        end
        
        % Make sure all times are relative to same reference point
        starttime = S1Ztstart;
        S1Zt = S1Zt + seconds(S1Ztstart-starttime);


    % Check to make sure there's actual data
    zeroind1 = find(S1Zraw == 0);
    if length(zeroind1) == length(S1Zraw)
        disp('All zeros!');
        continue
    end

    if(length(S1Zt)==0)
        display(['no data for ! station ',sta1]);
        continue
    end

    % Determine the time span to cut to ... this will change with
    % different segments
    clear tcut
    minT1Z = min(S1Zt);

    if length(S1Zraw) < 30000 
        disp(['Sta1 ',sta1,' : ',num2str(length(S1Zraw)),' is too short!'])
        continue
    end


        if(~exist('lat1','var'));

            lat1=S1.STLA;
            lon1=S1.STLO;
            dep1=S1.STEL; % depth is negative for OBS and positive for land stations

        end % if lat variabls

        stapairsinfo.stanames = {sta1};
        stapairsinfo.lats = [lat1];
        stapairsinfo.lons = [lon1];
        stapairsinfo.dt = dt;
        
%             % Frequency-time normalization
%             if IsFTN
%                 [ S1Zraw ] = FTN( S1Zraw, frange_FTN, dt );
%                 [ S2Zraw ] = FTN( S2Zraw, frange_FTN, dt );
%                 [ S1H1raw ] = FTN( S1H1raw, frange_FTN, dt );
%                 [ S2H1raw ] = FTN( S2H1raw, frange_FTN, dt );
%                 [ S1H2raw ] = FTN( S1H2raw, frange_FTN, dt );
%                 [ S2H2raw ] = FTN( S2H2raw, frange_FTN, dt );
%             end
        

        % START WINDOWING
        hour_length = winlength;

        nwin = floor(24/hour_length)*2-1; %
        win_length = hour_length*60*60/dt; % length of individual windows.
        win_start = 1;
        last_pt = win_length*.5*(nwin-1)+1+Nstart+win_length;
        if last_pt < length(S1Zraw)
            nwin = nwin + 1;
        end
		
%             tic
%             parfor iwin = 1:nwin
        ii_fftwin = 0;
        for iwin = 1:nwin
%				clear tcut S1R S2R S1T S2T S1Z S2Z fftS1R fftS2R fftS1T fftS2T fftS1Z fftS2Z
            clear fftS1Z_win_sum

			% cut in time
            if hour_length == 24
                pts_begin = Nstart;
                pts_end = length(S1Zraw)-Nstart;
            else
                pts_begin = win_length*.5*(iwin-1)+1+Nstart;
                pts_end = pts_begin+win_length;
            end

            if pts_begin > length(S1Zraw) || pts_end > length(S1Zraw)
				% disp('(Z) Points greater than the data... fixing window');
				pts_begin = length(S1Zraw)-win_length-Nstart;
                pts_end = pts_begin+win_length;
                %continue
            end
            tcut = [pts_begin:pts_end] * dt;

            % cut in tim Z for STA1
            S1Z=interp1(S1Zt,S1Zraw,tcut);
            S1Z(isnan(S1Z))=0;

            %detrend
            if IsDetrend
                S1Z = detrend(S1Z);
            end

            % Apply cosine taper
            if IsTaper
                S1Z = cos_taper(S1Z);
            end
            
            % Apply prefilter
            if IsPrefilter
                [bb,aa] = butter(2,frange_prefilt*2*dt);
                S1Z =  FiltFiltM(bb,aa,S1Z);
            end

            fftS1Z = fft(S1Z);
            if(~exist('fftS1Z_win_sum','var'))
                fftS1Z_win_sum = zeros(size(fftS1Z));
            end
            fftS1Z_win_sum = fftS1Z_win_sum + fftS1Z;
            ii_fftwin = ii_fftwin + 1;

        end % end window
        fftS1Z_day = fftS1Z_win_sum / ii_fftwin;
%             toc
        
        if(~exist('fftS1Z_day_sum','var'))
            fftS1Z_day_sum = zeros(size(fftS1Z_day));
        end
        fftS1Z_day_sum = fftS1Z_day_sum + fftS1Z_day;            
        ii_fftday = ii_fftday + 1;
        
        if(~exist('fftS1Z_month_sum','var'))
            fftS1Z_month_sum = zeros(size(fftS1Z_day));
        end
        fftS1Z_month_sum = fftS1Z_month_sum + fftS1Z_day;
        ii_fftday_month = ii_fftday_month + 1;
        
        if IsOutputDaystack
            % Save day stack
            daystr = datestr(starttime,'YYYYmmddHHMMSS');
            psdZ_daystack_path = [psd_singlestack_path,'psd',strNAMEcomp,'/'];
            clear fftS
            fftS = fftS1Z_day;
            save(sprintf('%s%s/%s_%s_f.mat',psdZ_daystack_path,sta1,sta1,daystr),'fftS','ii_fftwin','stapairsinfo','starttime');
        end
        if IsOutputMonthstack
            % Save 30 day (month) stack
            if month_counter == 30
                imonth = imonth + 1;
                psdZ_monthstack_path = [psd_monthstack_path,'psd',strNAMEcomp,'/'];
                clear fftS
                fftS1Z_month = fftS1Z_month_sum / ii_fftday_month;
                fftS = fftS1Z_month;
                save(sprintf('%s%s/%s_month%d_f.mat',psdZ_monthstack_path,sta1,sta1,imonth),'fftS','ii_fftday_month','stapairsinfo');
                month_counter = 0; % start over
                ii_fftday_month = 0;
            end
        end
    end % end hday
    
    if ii_fftday == 0
        display(['No data for ',sta1]);
        continue
    end
    fftS1Z_full = fftS1Z_day_sum / ii_fftday;

    if ii_fftwin > 1
        if IsFigure1
            f101 = figure(101);clf;
%                 set(gcf,'position',[400 400 600 300]);

            subplot(1,1,1);
            box on; hold on;
            T = length(fftS1Z_full);
            faxis = [0:(T-mod(T-1,2))/2 , -(T-mod(T,2))/2:-1]/dt/T;
            ind = find(faxis>0);
            plot(faxis(ind),abs(fftS1Z_full(ind)),'-','color',[0.7 0.7 0.7],'linewidth',2);
            plot(faxis(ind),smooth(abs(fftS1Z_full(ind)),100),'-k','linewidth',2);
            title(sprintf('%s PSD %s',sta1,strNAMEcomp(1)));
%                 xlim([0.01 1/(dt*2.5)])
            xlim([0.01 0.5])
            %xlim([0.04 0.16])
            xlabel('Frequency')
            set(gca,'xscale','log','YScale','log','fontsize',15,'linewidth',1.5)

%             print(f101,'-dpsc',[fig_winlength_path,sta1,'_',strNAMEcomp,'.ps']);
            saveas(f101, [fig_winlength_path,sta1,'_',strNAMEcomp,'.png']);
            %pause;
        end
        if IsOutputFullstack
            psdZ_fullstack_path = [psd_fullstack_path,'psd',strNAMEcomp,'/'];
            clear fftS
            fftS = fftS1Z_full;
            save(sprintf('%s%s/%s_f.mat',psdZ_fullstack_path,sta1,sta1),'fftS','ii_fftday','stapairsinfo');
        end
    end
end % ista1
