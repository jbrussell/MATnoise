% Convert matnoise data structure to the earthquake data structures used in
% ASWMS. Specifically, phase delay times determined from Bessel fitting of 
% ambient noise will be treated equivalently to the phase delay times 
% measured from GSDF. In other words, we assume that each station is a 
% virtual "earthquake" recorded by all the other stations.
%
clear

setup_parameters
setup_ErrorCode

compnoise = parameters.compnoise;
xspdir = parameters.xspdir; % 'phv_dir'; %'Nomelt3inttaper_iso.s0to333_br1avg'; %'4.0_S1_10pers_avg'; %'Nomelt3inttaper_iso.s0to333_br1avg'; %'4.0_S0_waverage';
windir = parameters.windir; %'window3hr'; 
N_wl = parameters.N_wl;
frange = parameters.frange; %[1/10 1/5]; % [Hz]
per_ind = parameters.per_ind; % [1:12]; % index of periods to consider

comp = parameters.component;

workingdir = parameters.workingdir;
CSoutputpath = [workingdir,'CSmeasure/'];

if ~exist(CSoutputpath,'dir')
	mkdir(CSoutputpath)
end

% Load station info
station_list = parameters.station_list;
[sta.name, sta.lat, sta.lon, sta.dep] = textread(station_list,'%s %f %f %f');

Xsp_path = ['../Xsp/',windir,'/fullStack/Xsp',compnoise{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',num2str(N_wl),'wl_',xspdir,'/'];
xspfiles = dir([Xsp_path,'*_xsp.mat']);

if ~exist(workingdir)
    mkdir(workingdir);
end

%% Loop over all xsp files and get virtual source and station names
str_src = {};
str_rec = {};
for ixsp = 1:length(xspfiles)
    toks = strsplit(xspfiles(ixsp).name,'_');
    str_src{ixsp} = toks{1};
    str_rec{ixsp} = toks{2};
end

%%
% Loop over virtual sources
periods_all = [];
for iev = 1:length(sta.name)
    event.id = sta.name{iev};
    stasrc = event.id;
    disp(['Working on virtual source: ',event.id]);
    
%     % Index the virtual source
%     Isrc = find(strcmp(str_src,event.id));
%     if isempty(Isrc)
%         disp(['No virtual source found: ',event.id])
%         continue
%     end
    
    % Loop through virtual source and
    % Form differential phase delays between pairs of stations
    stlas=[]; stlos=[];
    stnms={}; avgphv=[]; dists=[];
    CS = [];
    ii = 0;
    for ista1 = 1:length(sta.name)
        sta1 = sta.name{ista1};
        if strcmp(stasrc,sta1)
            continue
        end
        temp = load([Xsp_path,'/',stasrc,'_',sta1,'_xsp.mat']);
        xspinfo1 = temp.xspinfo;
        % Treat station 1 as the source
        evla = xspinfo1.lat1;
        evlo = xspinfo1.lon1;
        % Treat station 2 as the receiver (1)
        lat1 = xspinfo1.lat2;
        lon1 = xspinfo1.lon2;
        dist1 = vdist(evla,evlo,lat1,lon1)/1e3;
        per_ind1 = find(ismember(1:length(xspinfo1.per),per_ind));
        t1 = xspinfo1.tw(per_ind1);
        periods1 = xspinfo1.per(per_ind1);
        
        for ista2 = 1:length(sta.name)
            sta2 = sta.name{ista2};
            if strcmp(sta1,sta2) || strcmp(stasrc,sta2)
                continue
            end
            ii = ii+1;
            temp = load([Xsp_path,'/',stasrc,'_',sta2,'_xsp.mat']);
            xspinfo2 = temp.xspinfo;
            % Treat station 2 as the receiver (2)
            lat2 = xspinfo2.lat2;
            lon2 = xspinfo2.lon2;
            dist2 = vdist(evla,evlo,lat2,lon2)/1e3;
            per_ind2 = find(ismember(1:length(xspinfo2.per),per_ind));
            t2 = xspinfo2.tw(per_ind2);
            periods2 = xspinfo2.per(per_ind2);
            
            % Make sure same periods are indexed for both sta1 and sta2
            [per_ind12,I1,I2] = intersect(per_ind1,per_ind2);
            t1_persame = t1(I1);
            t2 = t2(I2);
            periods = periods2(I2);
            
            CS(ii).sta1 = find(strcmp(sta.name,sta1));
            CS(ii).sta2 = find(strcmp(sta.name,sta2));
%             CS(ii).win_cent_t = [];
            CS(ii).ddist = dist1-dist2;
            CS(ii).dtp = zeros(1,length(per_ind));
            CS(ii).amp = ones(1,length(per_ind));
            CS(ii).w = zeros(1,length(per_ind));
            CS(ii).exitflag = zeros(1,length(per_ind));
            CS(ii).isgood = zeros(1,length(per_ind));
            CS(ii).fiterr = zeros(1,length(per_ind));
            CS(ii).cohere = ones(1,length(per_ind));
            
            CS(ii).dtp(I2) = t1_persame-t2;
            CS(ii).amp(I2) = 1;
            CS(ii).w(I2) = 2*pi./periods;
            CS(ii).exitflag(I2) = 3;
            CS(ii).isgood(I2) = 1;
            
            periods_all = unique([periods_all(:)' periods(:)']);
        end
        % Source receiver distance for each station
        stlas(ista1) = xspinfo1.lat2;
        stlos(ista1) = xspinfo1.lon2;
        dists(ista1) = vdist(evla,evlo,stlas(ista1),stlos(ista1))/1e3;
    end
    % Treat sta1 as virtual source
    event.evla = evla;
    event.evlo = evlo;
    event.evdp = 0;
    event.Mw = 9999;
    
    % Get average phase velocity across the array and remove the outliers.
	clear avgphv
	for ip=1:length(periods)
		clear ddist dtp isgood
		for ics = 1:length(CS)
			ddist(ics) = CS(ics).ddist;
			dtp(ics) = CS(ics).dtp(ip);
			isgood(ics) = CS(ics).isgood(ip);
		end % end of ics
		goodind = find(isgood > 0);
		para = polyfit(ddist(goodind),dtp(goodind),1);
		err = abs(ddist*para(1) + para(2) - dtp);
		for ics = 1:length(CS)
			if err(ics) > parameters.tp_tol && CS(ics).isgood(ip) > 0
				CS(ics).isgood(ip) = ErrorCode.high_tp_err;
			end
			isgood(ics) = CS(ics).isgood(ip);
		end
		goodind = find(isgood > 0);
		para = polyfit(ddist(goodind),dtp(goodind),1);
		avgphv(ip) = 1./para(1);
	end % end of periods
    
    eventcs.CS = CS;
%     eventcs.autocor = event.autocor;
    eventcs.id = event.id;
    eventcs.avgphv = avgphv;
    eventcs.stlas = stlas;
    eventcs.stlos = stlos;
    eventcs.stnms = stnms;
    eventcs.evla = event.evla;
    eventcs.evlo = event.evlo;
    eventcs.evdp = event.evdp;
    eventcs.dists = dists; % km
    eventcs.eventmatfile = ''; %[eventmatpath,matfiles(ie).name];
    eventcs.Mw = event.Mw;

    matfilename = [CSoutputpath,char(event.id),'_cs_',comp,'.mat'];
	save(matfilename,'eventcs')
	disp(['Save to ',matfilename]);
end

periods = periods_all;
save('periods.mat','periods');