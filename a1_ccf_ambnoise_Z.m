% Calculate ambient noise cross correlation record from multiple stationpairs 
% for Vertical (Z) only using the methods from Bensen et al. (2007) GJI 
% DOI:10.1111/j.1365-246X.2007.03374.x
%
% Expects files organized like so:
% {datadirectory}/{station}/{station}.{yyyy}.{jday}.{hh}.{mm}.{SS}.{COMP}.sac
%  e.g.: mydata/CC05/CC05.2018.112.00.00.00.BDH.sac
%
% JBR, Jan 2020: Implemented frequency-time normalization after 
% Shen et al. (2012) BSSA; DOI:10.1785/0120120023. This greatly improves signal
% extraction compared to typical one-bit noralization and whitening of Bensen et
% al. (2007) GJI. Faster FiltFiltM() can be replaced with MATLAB's slower 
% built-in filtfilt().
%
% JBR, update: We have found that doing no time or frequency normalization at all
% can produce higher SNR traces than doing one-bit or time-frequency normalization. 
% Therefore, the default is to use the raw seismograms as is without any preprocessing.
%
% (NOTE: FUNCTIONIZE IN THE FUTURE)
% Patty Lin -- 10/2014
% Natalie Accardo
% Josh Russell
% https://github.com/jbrussell
clear;
setup_parameters;

strSACcomp = 'Z';
strNAMEcomp = 'ZZ';
IsFigure1 = 1;
IsFigure2 = 0;

% OUTPUT SETTINGS
IsOutputFullstack = 1; % Save full year ccf stacks
IsOutputMonthstack = 0; % save month ccf stacks
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
figpath = parameters.figpath;
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
ccf_path = parameters.ccfpath;
ccf_winlength_path = [ccf_path,'window',num2str(winlength),'hr/'];
ccf_singlestack_path = [ccf_winlength_path,'single/'];
ccf_daystack_path = [ccf_winlength_path,'dayStack/'];
ccf_monthstack_path = [ccf_winlength_path,'monthStack/'];
ccf_fullstack_path = [ccf_winlength_path,'fullStack/'];

if ~exist(ccf_path)
    mkdir(ccf_path)
end
if ~exist(ccf_winlength_path)
    mkdir(ccf_winlength_path)
end
if ~exist(ccf_singlestack_path)
    mkdir(ccf_singlestack_path)
end
if ~exist(ccf_daystack_path)
    mkdir(ccf_daystack_path)
end
if ~exist(ccf_monthstack_path)
    mkdir(ccf_monthstack_path)
end
if ~exist(ccf_fullstack_path)
    mkdir(ccf_fullstack_path)
end

PATHS = {ccf_singlestack_path; ccf_daystack_path; ccf_monthstack_path; ccf_fullstack_path};
for ipath = 1:length(PATHS)
    ccfZ_path = [PATHS{ipath},'ccf',strNAMEcomp,'/'];
    if ~exist(ccfZ_path)
        mkdir(ccfZ_path);
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
nsta=parameters.nsta; % number of target stations to calculate for

% READ OBS ORIENTATIONS
[slist, orientations] = textread(orientation_path,'%s%f\n');

% Calculate filter coefficients for FTN
if IsFTN
    [ b, a ] = get_filter_TFcoeffs( frange_FTN, dt );
end

for ista1=1:nsta

    sta1=char(stalist(ista1,:));
    % Build station directories
    for ipath = 1:length(PATHS)
        ccfZ_path = [PATHS{ipath},'ccf',strNAMEcomp,'/'];
        if ~exist([ccfZ_path,sta1])
            mkdir([ccfZ_path,sta1]);
        end
    end
    seisZ_path = [seis_winlength_path,strNAMEcomp(1),'/'];
    if ~exist([seisZ_path,sta1])
        mkdir([seisZ_path,sta1]);
    end

    list1 = dir([datadir,sta1,'/*',strSACcomp,'.sac']);

    for ista2=1:nsta
        clear lat1 lat2 lon1 lon2 dist az baz vec_tz2 Z2raw vec_tz Z1raw

        sta2=char(stalist(ista2,:));

        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end


        % check to see if we've already done this ccf
        if exist([ccfZ_path,sta1,'/',sta1,'_',sta2,'_f.mat'])
            display('CCF already exist, skip this pair');
            continue
        end

        display(['performing cross-correlation for staion pair : ',sta1,'  ', sta2]);
        % -------------loop through each half day--------------------
        nday_stack=0;
        coh_sumZ = 0;
        coh_num = 0;

        % Get a list of all available data
        ihday = 0;
        month_counter = 0;
        imonth = 0;
        for ifil = 1:length(list1)
            file1cZ = list1(ifil).name;

            % Check that day file exists for station 2
            Nchar = length(sta1);
            file2cZ = dir([datadir,sta2,'/',sta2,file1cZ(Nchar+1:end)]);
            str = strsplit(file1cZ,'.');
            hdayid = [str{2},'.',str{3},'.',str{4},'.',str{5},'.',str{6}];
            
            if isempty(file2cZ)
                disp(['No data for ',sta2,' on day ',hdayid,'... skipping'])
                continue
            end
            file2cZ = file2cZ.name;

            if month_counter == 0
                coh_sumZ_month = 0;
                coh_num_month = 0;
            end
            clear data1cZ data2cZ
            ihday = ihday +1;
            month_counter = month_counter + 1;
            clear temp
            %temp = strsplit(daylist1(ihday).name,'.');

            disp(['Looking at ',hdayid,' ',sta2]);

            data1cZ= dir([datadir,sta1,'/',year,'/',sta1,'.',hdayid,'.*',strSACcomp,'.sac']);
            data2cZ= dir([datadir,sta2,'/',year,'/',sta2,'.',hdayid,'.*',strSACcomp,'.sac']);

            data1cZ =  [datadir,sta1,'/',year,'/',data1cZ.name];
            data2cZ =  [datadir,sta2,'/',year,'/',data2cZ.name];

            %------------------- TEST IF DATA EXIST------------------------
            [S1Zt,S1Zraw,S1,S1Ztstart] = load_sac(data1cZ);
            [S2Zt,S2Zraw,S2,S2Ztstart] = load_sac(data2cZ);
            
            % Check that sample rates are the same
            if S1.DELTA ~= S2.DELTA
                error('S1 and S2 sample rates don''t match!');
            end
            
            % Make sure all times are relative to same reference point
            starttime = S1Ztstart;
            S1Zt = S1Zt + seconds(S1Ztstart-starttime);
            S2Zt = S2Zt + seconds(S2Ztstart-starttime);
            
            % Ensure that files have same start time to within 1 sample
            if abs(seconds(S1Ztstart-S2Ztstart)) > S1.DELTA
                error('Station files do not have same start time');
            end
            
            % Make sure sample rates all match
            if (abs(S1.DELTA-dt) >= 0.01*dt ) || (abs(S2.DELTA-dt) >= 0.01*dt )
                error('sampling interval does not match data! check dt');
            end

            %------------------- Remove instrument response ------------------------
        if IsRemoveIR
            pzfile1 = dir([PZpath,'/RESP.*.',sta1,'.*.*Z']); % PZ for H1 and H2 are identical
            pzfile2 = dir([PZpath,'/RESP.*.',sta2,'.*.*Z']);

            % Read RESP file for station 1
            [z,p,c,units] = read_sac_RESP([PZpath,pzfile1.name],units_RemoveIR);

            dt1 = abs(S1Zt(1)-S1Zt(2));
            dt2 = abs(S2Zt(1)-S2Zt(2));

            % Remove instrument response for station 1 Z
            S1Zraw = rm_resp(S1Zraw,z,p,c,dt1);

            % Read RESP file for station 2
            [z,p,c,units] = read_sac_RESP([PZpath,pzfile2.name],units_RemoveIR);

            % Remove instrument response for station 2 Z
            S2Zraw = rm_resp(S2Zraw,z,p,c,dt2);
        end


        % Check to make sure there's actual data
        zeroind1 = find(S1Zraw == 0);
        zeroind2 = find(S2Zraw == 0);
        if length(zeroind1) == length(S1Zraw) || length(zeroind2) == length(S2Zraw)
            disp('All zeros!');
            continue
        end

        if(length(S1Zt)*length(S2Zt)==0)
            display(['no data for ! station ',sta2]);
            continue
        end

        % Determine the time span to cut to ... this will change with
        % different segments
        clear tcut
        minT1Z = min(S1Zt);
        minT2Z = min(S2Zt);

        if length(S1Zraw) < 30000 
            disp(['Sta1 ',sta1,' : ',num2str(length(S1Zraw)),' is too short!'])
            continue
        elseif length(S2Zraw) < 30000 
            disp(['Sta2 ',sta2,' : ',num2str(length(S2Zraw)),' is too short!'])
            continue
        end


            if(~exist('lat2','var'));

                lat1=S1.STLA;
                lon1=S1.STLO;
                dep1=S1.STEL; % depth is negative for OBS and positive for land stations


                lat2=S2.STLA;
                lon2=S2.STLO;
                dep2=S2.STEL; % depth is negative for OBS and positive for land stations


                % Get the interstation distance and azimuth
                [delta,S1az]=distance(lat1,lon1,lat2,lon2);
                [delta,S2az]=distance(lat2,lon2,lat1,lon1);

                dist=deg2km(delta);

                if(dist < dist_min)
                    display(['distance shorter than ',num2str(dist_min),' km, skip']);
                    break
                end
            end % if lat variabls

            stapairsinfo.stanames = {sta1,sta2};
            stapairsinfo.lats = [lat1,lat2];
            stapairsinfo.lons = [lon1,lon2];
            stapairsinfo.dt = dt;
            stapairsinfo.r = dist;
            
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
            coh_sumZ_day = 0;
            coh_num_day = 0;
            last_pt = win_length*.5*(nwin-1)+1+Nstart+win_length;
            if last_pt < length(S1Zraw)
                nwin = nwin + 1;
            end
			
%             tic
            parfor iwin = 1:nwin
%				clear tcut S1R S2R S1T S2T S1Z S2Z fftS1R fftS2R fftS1T fftS2T fftS1Z fftS2Z

				% cut in time
                if hour_length == 24
                    pts_begin = Nstart;
                    pts_end = length(S1Zraw)-Nstart;
                else
                    pts_begin = win_length*.5*(iwin-1)+1+Nstart;
                    pts_end = pts_begin+win_length;
                end

                if pts_begin > length(S1Zraw) || pts_begin > length(S2Zraw) || pts_end > length(S1Zraw) || pts_end > length(S2Zraw)
					% disp('(Z) Points greater than the data... fixing window');
					pts_begin = length(S1Zraw)-win_length-Nstart;
                    pts_end = pts_begin+win_length;
                    %continue
                end
                tcut = [pts_begin:pts_end] * dt;

                % cut in tim Z for STA1
                S1Z=interp1(S1Zt,S1Zraw,tcut);
                S1Z(isnan(S1Z))=0;

                % cut in tim Z for STA2
                S2Z=interp1(S2Zt,S2Zraw,tcut);
                S2Z(isnan(S2Z))=0;

                %detrend
            if IsDetrend
                S1Z = detrend(S1Z);
                S2Z = detrend(S2Z);
            end

            % Apply cosine taper
            if IsTaper
                S1Z = cos_taper(S1Z);
                S2Z = cos_taper(S2Z);
            end
            
            % Apply prefilter
            if IsPrefilter
                [b,a] = butter(2,frange_prefilt*2*dt);
                S1Z =  FiltFiltM(b,a,S1Z);
                S2Z =  FiltFiltM(b,a,S2Z);
            end

                if IsFigure2
                    figure(49)
                    clf

                    %Z
                    subplot(5,1,1)
                    plot(tcut,S1Z,'-k')
                    %ylim([-0.15e-5 0.15e-5])
                    xlim([tcut(1) tcut(end)])
                    title(strNAMEcomp(1));
                    hold on

                    pause;
                    %return
                end

                %-------------------- Vertical Component --------------
                if IsFTN
                    % Frequency-time normalization
                    [ S1Z ] = FTN( S1Z, b, a  );
                    [ S2Z ] = FTN( S2Z, b, a  );
                    fftS1Z = fft(S1Z);
                    fftS2Z = fft(S2Z);
                else
                    % One-bit normalization
                    if IsOBN
                        S1Z = runwin_norm(S1Z);
                        S2Z = runwin_norm(S2Z);
                    end
                    %fft
                    fftS1Z = fft(S1Z);
                    fftS2Z = fft(S2Z);
                    %Whiten
                    if IsSpecWhiten
                        fftS1Z = spectrumwhiten_smooth(fftS1Z,0.001);
                        fftS2Z = spectrumwhiten_smooth(fftS2Z,0.001);
                    end
                end

                % calcaulate daily CCF and stack for radial
                coh_trace = fftS1Z .* conj(fftS2Z);
                coh_trace = coh_trace ./ abs(fftS1Z) ./ abs(fftS2Z);
                coh_trace(isnan(coh_trace)) = 0;
                coh_sumZ = coh_sumZ + coh_trace;
                coh_trace_Z = coh_trace;
                coh_sumZ_day = coh_sumZ_day + coh_trace;
                coh_sumZ_month = coh_sumZ_month + coh_trace;

                % coh_num = coh_num + 1;
                coh_num_day = coh_num_day + 1;
                coh_num_month = coh_num_month + 1;
    %             toc

                if IsOutputSinglestack % save individual xcor
                    ccfZ_singlestack_path = [ccf_singlestack_path,'ccf',strNAMEcomp,'/'];
                    save(sprintf('%s%s/%s_%s_%d_f.mat',ccfZ_singlestack_path,sta1,sta1,sta2,coh_num),'coh_trace_Z','stapairsinfo');
                end
                if IsOutputSeismograms % save seismograms before xcor
                    seisZ_path = [seis_winlength_path,strNAMEcomp(1),'/'];
                    save(sprintf('%s%s/%s_%d_f.mat',seisZ_path,sta1,sta1,coh_num),'S1Z','stapairsinfo');
                end
            end % end window
%             toc
            coh_num = coh_num + coh_num_day;
            
            if IsOutputDaystack
                % Save day stack
                daystr = datestr(starttime,'YYYYmmddHHMMSS');
                ccfZ_daystack_path = [ccf_daystack_path,'ccf',strNAMEcomp,'/'];
                clear coh_sum
                coh_sum = coh_sumZ_day;
                save(sprintf('%s%s/%s_%s_%s_f.mat',ccfZ_daystack_path,sta1,sta1,sta2,daystr),'coh_sum','coh_num_day','stapairsinfo','starttime');
            end
            if IsOutputMonthstack
                % Save 30 day (month) stack
                if month_counter == 30
                    imonth = imonth + 1;
                    ccfZ_monthstack_path = [ccf_monthstack_path,'ccf',strNAMEcomp,'/'];
                    clear coh_sum
                    coh_sum = coh_sumZ_month;
                    save(sprintf('%s%s/%s_%s_month%d_f.mat',ccfZ_monthstack_path,sta1,sta1,sta2,imonth),'coh_sum','coh_num_month','stapairsinfo');
                    month_counter = 0; % start over
                end
            end
        end % end hday

        if coh_num > 1
            if IsFigure1
                f101 = figure(101);clf;
%                 set(gcf,'position',[400 400 600 300]);

                subplot(3,1,3)
                T = length(coh_sumZ);
                faxis = [0:(T-mod(T-1,2))/2 , -(T-mod(T,2))/2:-1]/dt/T;
                ind = find(faxis>0);
                plot(faxis(ind),smooth(real(coh_sumZ(ind)/coh_num),100));
                title(sprintf('%s %s coherency %s ,station distance: %f km',sta1,sta2,strNAMEcomp(1),dist));
                xlim([0.01 0.5])
                %xlim([0.04 0.16])
                xlabel('Frequency')
                drawnow

                print(f101,'-dpsc',[fig_winlength_path,sta1,'_',sta2,'_',strNAMEcomp,'.ps']);
                %pause;
            end
            if IsOutputFullstack
                ccfZ_fullstack_path = [ccf_fullstack_path,'ccf',strNAMEcomp,'/'];
                clear coh_sum
                coh_sum = coh_sumZ;
                save(sprintf('%s%s/%s_%s_f.mat',ccfZ_fullstack_path,sta1,sta1,sta2),'coh_sum','coh_num','stapairsinfo');
            end
        end
    end % ista2

end % ista1

delete(gcp('nocreate')); % remove parallel pools