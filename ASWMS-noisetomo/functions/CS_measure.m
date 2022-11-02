function CS = CS_measure(event,sta1,sta2,parameters)
% Main function to perform GSDF measurement
%
	setup_ErrorCode;
	isdebug = parameters.isdebug;
    isfigure = parameters.isfigure;
	is_winx2 = parameters.is_winx2;

	refv = parameters.refv;
	refphv = parameters.refphv;
	periods = parameters.periods;
	min_width = parameters.min_width;
	max_width = parameters.max_width;
	wintaperlength = parameters.wintaperlength;
	prefilter = parameters.prefilter;
	xcor_win_halflength = parameters.xcor_win_halflength;
	Nfit = parameters.Nfit;
	Ncircle = parameters.Ncircle;
	xcor_win_iter = parameters.xcor_win_iter;
	if length(xcor_win_iter) ~= length(periods)
		disp('xcor_win_iter should be the same size as periods!');
		if length(xcor_win_iter) < length(periods)
			xcor_win_iter = [zeros(1,length(periods) - length(xcor_win_iter)),xcor_win_iter];
		end
	end

	CS = init_CSstruct;
	v1 = event.winpara(1); t1=event.winpara(2); v2=event.winpara(3); t2=event.winpara(4);

	CS.sta1 = sta1;
	CS.sta2 = sta2;
    if isdebug
        disp([num2str(sta1),num2str(sta2)]);
    end
	if length(event.stadata(sta1).data) < 20 || length(event.stadata(sta2).data) < 20
		CS.exitflag(:) = ErrorCode.sta_lackdata;
		return
	end

	% read in data for station 1 and apply prefilter
	data1 = event.stadata(sta1).data;
	bgtime = event.stadata(sta1).otime - event.otime;
	dt1 = event.stadata(sta1).delta;
	Nt = length(event.stadata(sta1).data);
	fN = 1/2/dt1;
	[b,a] = butter(2,[1/prefilter(2)/fN, 1/prefilter(1)/fN]);
	data1 = filtfilt(b,a,data1);
	taxis1 = bgtime + [0:Nt-1]'*dt1;
	dist1 = event.stadata(sta1).dist;
	winbgt = dist1/v1+t1;
	winendt = dist1/v2+t2;
	if taxis1(1) > winbgt || taxis1(end) < winendt 
		CS.exitflag(:) = ErrorCode.sta_lackdata;
		return
	end

	% read in data for station 2 and apply prefilter
	data2 = event.stadata(sta2).data;
	bgtime = event.stadata(sta2).otime - event.otime;
	dt2 = event.stadata(sta2).delta;
	Nt = length(event.stadata(sta2).data);
	fN = 1/2/dt2;
	[b,a] = butter(2,[1/prefilter(2)/fN, 1/prefilter(1)/fN]);
	data2 = filtfilt(b,a,data2);
	taxis2 = bgtime + [0:Nt-1]'*dt2;
	dist2 = event.stadata(sta2).dist;
	winbgt = dist2/v1+t1;
	winendt = dist2/v2+t2;
	if taxis2(1) > winbgt || taxis2(end) < winendt 
		CS.exitflag(:) = ErrorCode.sta_lackdata;
		return
	end

	% resample the data if necessary
	if dt1 > dt2
		new_taxis2 = taxis2(1):dt1:taxis2(end);
		data2 = interp1(taxis2,data2,new_taxis2);
		taxis2 = new_taxis2;
		dt2 = dt1;
	elseif dt1 < dt2
		new_taxis1 = taxis1(1):dt2:taxis1(end);
		data1 = interp1(taxis1,data1,new_taxis1);
		taxis1 = new_taxis1;
		dt1 = dt2;
	end

	% window data2
	winbgt = dist2/v1+t1;
	winendt = dist2/v2+t2;
	win_data2 = flat_hanning_win(taxis2,data2,winbgt,winendt,wintaperlength);
	
	if is_winx2
		% window data1
		winbgt = dist1/v1+t1;
		winendt = dist1/v2+t2;
		data1 = flat_hanning_win(taxis1,data1,winbgt,winendt,wintaperlength);
	end

	% apply cross-correlation
	[xcor,lag] = xcorr(data1,win_data2,...
        floor(10*max(periods)/dt1+abs(taxis1(1)-taxis2(1))));
	lag = lag.*dt1;
	lag = lag + taxis1(1) - taxis2(1);

	if isfigure
		figure(43)
		clf
		subplot(3,1,1)
		plot(taxis1,data1);
		xlim([0 dist2/2])
		subplot(3,1,2)
		plot(taxis2,win_data2);
		xlim([0 dist2/2])
		subplot(3,1,3)
		plot(lag,xcor);
		xlim([-1000 1000])
	end

	%Find the window center (max amplitude within the window)
	win_cent_t = (dist1-dist2)/refv;
	search_win_ind = find( lag > win_cent_t-xcor_win_halflength &...
		lag < win_cent_t + xcor_win_halflength );
	[max_xcor_amp win_cent_i] = max(xcor(search_win_ind));
	win_cent_i = search_win_ind(win_cent_i);
	win_cent_t = lag(win_cent_i);
	CS.win_cent_t = win_cent_t;
	CS.ddist = dist1 - dist2;
	% apply the window function
	%win_xcor = hanning_win(lag,xcor,win_cent_t,xcor_win_halflength*2);
	win_xcor = flat_hanning_win(lag,xcor,win_cent_t-xcor_win_halflength,win_cent_t+xcor_win_halflength,round(xcor_win_halflength/2));


	if isfigure
		figure(44)
		clf
		subplot(2,1,1)
		plot(lag,xcor);
		xlim([-500 500])
		subplot(2,1,2)
		plot(lag,win_xcor);
		xlim([-500 500])
	end

	% Apply Narrow-band filter
	clear gaus_filters nband_win_xcors
	Nt = length(win_xcor);
	[gaus_filters,faxis] = build_gaus_filter(1./periods,dt1,Nt,min_width,max_width);
	fft_win_xcor = fft(win_xcor);
	if size(fft_win_xcor) == 1
		fft_win_xcor = fft_win_xcor';
	end
	for ip = 1:length(periods)
		nband = fft_win_xcor .* [gaus_filters(:,ip); zeros(Nt-length(gaus_filters(:,ip)),1)];
		nband = ifft(nband);
		nband = 2*real(nband);
		nband_win_xcors(:,ip) = nband;
	end % end of periods loop


	% fitting with five-parameter wavelet
	for ip = 1:length(periods)
		[para,resnorm,residual, exitflag] = gsdffit(nband_win_xcors(:,ip),lag,1./periods(ip),Nfit);
		CS.fitpara(:,ip) = para(:);
		CS.fiterr(ip) = resnorm./para(1)^2./Nfit./periods(ip);
		CS.dtp(ip) = para(4);
		CS.dtg(ip) = para(5);
		CS.amp(ip) = para(1);
		CS.w(ip) = para(2);
		CS.sigma(ip) = para(3);
		CS.exitflag(ip) = exitflag;
	end

	% Iteratively correct for windowing effect
	for ip = 1:length(periods)
		if xcor_win_iter(ip)
			% re-center the window
			win_cent_t = CS.dtg(ip);
			%win_xcor = hanning_win(lag,xcor,win_cent_t,xcor_win_halflength*2);
			win_xcor = flat_hanning_win(lag,xcor,win_cent_t-xcor_win_halflength,win_cent_t+xcor_win_halflength,round(xcor_win_halflength/2));

			fft_win_xcor = fft(win_xcor);
			if size(fft_win_xcor) == 1,fft_win_xcor = fft_win_xcor'; end
			% narrow-band filter
			nband = fft_win_xcor .* [gaus_filters(:,ip); zeros(Nt-length(gaus_filters(:,ip)),1)];
            nband = ifft(nband);
            nband = 2*real(nband);
            % fit the wavelet again
			[para,resnorm,residual, exitflag] = gsdffit(nband(:),lag,1./periods(ip),Nfit);
			CS.fitpara(:,ip) = para(:);
			CS.fiterr(ip) = resnorm./para(1)^2./Nfit./periods(ip);
			CS.dtp(ip) = para(4);
			CS.dtg(ip) = para(5);
			CS.amp(ip) = para(1);
			CS.w(ip) = para(2);
			CS.sigma(ip) = para(3);
			CS.exitflag(ip) = exitflag;
		end
	end

	% Correct for cycle skipping for dtp
	for ip = 1:length(periods)
		syndtp = CS.ddist./refphv(ip);
		testdtp = CS.dtp(ip) + [-Ncircle:Ncircle]*periods(ip);
		[temp besti] = min(abs(testdtp - syndtp));
		CS.dtp(ip) = testdtp(besti);
	end

	if isfigure
		figure(45)
		clf
		hold on
		[xi yi] = ndgrid(lag,periods);
		for ip = 1:length(periods)
			norm_nbands(:,ip) = nband_win_xcors(:,ip)./max(abs(nband_win_xcors(:,ip)));
		end
		contourf(xi,yi,norm_nbands);
		for ip=1:length(periods)
			plot(CS.dtp(ip),periods(ip),'kx','linewidth',2);
		end
		xlim([-3*max(periods) 3*max(periods)]);
%         	pause
        drawnow
    end
    
    if isfigure
        figure(46)
        clf;
        
        subplot(2,1,1);
        hold on;
        plot(periods,CS.fiterr,'o')
        xlabel('periods (s)');
        ylabel('fit err');
        
        if sta1 ~= sta2
            subplot(2,1,2);
            cohere = CS.amp.^2./event.autocor(sta1).amp./event.autocor(sta2).amp;
            hold on;
            plot(periods,cohere,'or');
            xlabel('periods (s)');
            ylabel('coherence');
        end
    end


end % end of function
