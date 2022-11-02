function outevent = exam_winpara(event)

setup_parameters;

lalim = parameters.lalim;
lolim = parameters.lolim;

outevent = event;

periods = parameters.periods;
stadata = event.stadata;
evla = event.evla;
evlo = event.evlo;

time_range = [0 6000];
N_trace = 20;
stlas = [stadata.stla];
stlos = [stadata.stlo];
[dists azi] = distance(evla,evlo,stlas,stlos);

corner_dists(1) = distance(evla,evlo,lalim(1),lolim(1));
corner_dists(2) = distance(evla,evlo,lalim(1),lolim(2));
corner_dists(3) = distance(evla,evlo,lalim(2),lolim(1));
corner_dists(4) = distance(evla,evlo,lalim(2),lolim(2));
dist_range = [min(corner_dists) max(corner_dists)];
%dist_range = [min(dists)-1 max(dists)+1];

% prepare the data
for ista=1:length(stadata);
	bgtime = stadata(ista).otime - event.otime;
	npts = length(stadata(ista).data);
	delta = stadata(ista).delta;
	timeaxis = bgtime:delta:bgtime+delta*(npts-1);
	stadata(ista).timeaxis = timeaxis;
	bpfilter = [periods(end) periods(1)];
	fN = 1./delta./2;
	w = 1./bpfilter./fN;
	[b a] = butter(2,w);
	stadata(ista).odata = stadata(ista).data;
	stadata(ista).data = filtfilt(b,a,stadata(ista).data);	
end
% parameters that not need to be changed.

ori_dist_range = dist_range;
ori_time_range = time_range;
azi_range = [min(azi) max(azi)];
zoom_level = 1;
hist_time_range(zoom_level,:) = time_range;
hist_dist_range(zoom_level,:) = dist_range;
freq_band = 0; 
comp = 1;
single_norm = 1;
amp = 5;
norm_amp = 1;
isfill = 0;
is_reduce_v = 0;
ref_v = 10;
is_dist = 1;
is_cheatsheet = 0;
is_bin = 1;
is_mark = 0;
amp_diff_tol = 5;


while 1


	figure(99)
	clf
	hold on
	set(gca,'YDir','reverse');
	clear max_amp
	for ista = 1:length(stadata)
		if dists(ista) < dist_range(1) || dists(ista) > dist_range(2)
			max_amp(ista) = NaN;
			continue;
		end
		timeaxis = stadata(ista).timeaxis;
		if is_reduce_v
			timeaxis = timeaxis - deg2km(dists(ista))./ref_v;
		end
		ind = find(timeaxis > time_range(1) & timeaxis < time_range(2));
		data = stadata(ista).data;
		data = data(ind);
		if ~isempty(data)
			max_amp(ista) = max(abs(data));
		else
			max_amp(ista) = NaN;
		end
	end
	norm_amp = nanmedian(max_amp);
	dist_bin = linspace(dist_range(1),dist_range(2),N_trace);
	plot_bin = zeros(size(dist_bin));
	ind = find(dists>dist_range(1) & dists < dist_range(2));
	azi_range = [min(azi(ind)) max(azi(ind))];
	azi_bin = linspace(azi_range(1),azi_range(2),N_trace);
	
	for ista = 1:length(stadata)
		if dists(ista) < dist_range(1) || dists(ista) > dist_range(2)
			continue;
		end
		timeaxis = stadata(ista).timeaxis;
		if is_reduce_v
			timeaxis = timeaxis - deg2km(dists(ista))./ref_v;
		end
		ind = find(timeaxis > time_range(1) & timeaxis < time_range(2));
		if isempty(ind)
			continue;
		end
		timeaxis = timeaxis(ind);
		data = stadata(ista).data;
		data = data(ind);
		if single_norm
			data = data./max(abs(data));
		else
			data = data./norm_amp;
			if max(abs(data)) > amp_diff_tol
				data(:) = 0;
			end
		end
		if is_dist
			if is_bin
				bin_id = round((dists(ista)-dist_bin(1))./(dist_bin(2)-dist_bin(1)));
				if bin_id == 0
					bin_id = bin_id+1;
				end
				if bin_id > length(plot_bin)
					bin_id = bin_id-1;
				end
				plot_bin(bin_id) = plot_bin(bin_id)+1;
				if plot_bin(bin_id) > 1
					continue;
				end
			end
			trace_amp = amp*diff(dist_range)/(2*N_trace);
			plot(timeaxis,data*trace_amp+dists(ista),'k');
			text(max(timeaxis),dists(ista),num2str(stadata(ista).snr));
			if isfill
				data(find(data > 0)) = 0;
				area(timeaxis,data*trace_amp+dists(ista),dists(ista),'facecolor','b');
			end
			if is_mark
				plot(markertime,markerdist,'m','linewidth',2);
			end
		else
			if is_bin
				bin_id = round((azi(ista)-azi_bin(1))./(azi_bin(2)-azi_bin(1)));
				if bin_id < 1 || bin_id > length(plot_bin)
					continue;
				end
				plot_bin(bin_id) = plot_bin(bin_id)+1;
				if plot_bin(bin_id) > 1
					continue;
				end
			end
			trace_amp = amp*diff(azi_range)/(2*N_trace);
			plot(timeaxis,data*trace_amp+azi(ista),'k');
			if isfill
				data(find(data > 0)) = 0;
				area(timeaxis,data*trace_amp+azi(ista),azi(ista),'facecolor','b');
			end
		end
	end % end of station loop
	win_T1 = deg2km(dist_range)./event.winpara(1)+event.winpara(2);
	win_T2 = deg2km(dist_range)./event.winpara(3)+event.winpara(4);
	plot(win_T1,dist_range,'r','linewidth',2);
	plot(win_T2,dist_range,'r','linewidth',2);
	if is_dist
		ylim(dist_range);
	else
		ylim(azi_range);
	end
	xlim(time_range);
	if is_reduce_v
		refvstr = [num2str(ref_v),' km/s'];
		raypstr = [num2str(1/km2deg(ref_v)), ' s/deg '];
	else
		refvstr = 'None';
		raypstr = 'None';
	end
	stemp = sprintf('Comp: %s, ampnorm: %d, isgood: %d, winpara: %5.1f %5.1f %5.1f %5.1f',parameters.component,single_norm,event.isgood,event.winpara(1),event.winpara(2),event.winpara(3),event.winpara(4));
	title(stemp,'fontsize',15);
	xlabel('Time /s','fontsize',15)
	if is_dist
		ylabel('Distance /degree','fontsize',15)
	else
		ylabel('Azimuth /degree','fontsize',15)
	end

	[x y bot] = ginput(1);
	if ~is_dist
		is_dist = 1;
		continue;
	end
	if bot == 'q'
		break;
	end
	if bot == 'f'
		[x2 y2 ip] = ginput(1);
		temp = ip-'0';
		if temp > 0 && temp <= length(periods);
			ip = temp;
			bpfilter = [periods(ip)*1.1 periods(ip)*0.9];
			fN = 1./delta;
			w = 1./bpfilter./fN;
			[b a] = butter(2,w);
			for ista=1:length(stadata)
				stadata(ista).data = filtfilt(b,a,stadata(ista).odata);	
			end
		end
		if temp == 0
			bpfilter = [periods(end) periods(1)];
			fN = 1./delta;
			w = 1./bpfilter./fN;
			[b a] = butter(2,w);
			for ista=1:length(stadata)
				stadata(ista).data = filtfilt(b,a,stadata(ista).odata);	
			end
		end
			
	end
	if bot == 'g'
		event.isgood = ~event.isgood;
	end
	if bot == 'd'
		is_dist = ~is_dist;
	end
	if bot == '.'
		if is_mark
			is_mark = 0;
		else
			is_mark = 1;
			[x2 y2] = ginput(1);
			markertime = [x x2];
			markerdist = [y y2];
		end
	end
	if bot == 'a'
		[x2 y2 ampstr] = ginput(1);
		temp = ampstr - '0';
		if temp > 0 & temp < 10
			amp = temp;
		end
	end
	if bot == 'i'
		rayp = input('Reduce slowness(s/deg):');
		new_ref_v = deg2km(1/rayp);
		if is_reduce_v
			time_range = time_range + deg2km(mean(dist_range))./ref_v - deg2km(mean(dist_range))./new_ref_v;;
			for izoom = 1:size(hist_time_range,1)
				hist_time_range(izoom,:) = hist_time_range(izoom,:) ...
					+ deg2km(mean(dist_range))./ref_v - deg2km(mean(dist_range))./new_ref_v;
			end
		end
		ref_v = new_ref_v;
	end
	if bot == 'n'
		single_norm = ~single_norm;
	end
	if bot == 'p'
		[x2 y2] = ginput(1);
		plot([x x2],[y y2],'b');
		times = [x x2];
		epi_dists = deg2km([y y2]);
		para = polyfit(epi_dists,times,1);
		para(1) = 1./para(1);
		para
		temp = input('T1 or T2 (1/2):');
		if temp == 1
			event.winpara(1) = para(1);
			event.winpara(2) = para(2);
		elseif temp == 2
			event.winpara(3) = para(1);
			event.winpara(4) = para(2);
		end
	end
	if bot == 'x'  % changing time range
		plot([x x],dist_range,'r','linewidth',2);
		[x2 y2] = ginput(1);
		plot([x2 x2],dist_range,'r','linewidth',2);
		pause(0.5);
		if x2 > x
			time_range = [x x2];
		else
			time_range = [x2 x];
		end
		zoom_level = zoom_level + 1;
		hist_time_range(zoom_level,:) = time_range;
		hist_dist_range(zoom_level,:) = dist_range;
	end
	if bot == 'y' % change distance range
		plot(time_range,[y y],'r','linewidth',2);
		[x2 y2] = ginput(1);
		plot(time_range,[y2 y2],'r','linewidth',2);
		pause(0.5)
		if y2 > y
			dist_range = [y y2];
		else
			dist_range = [y2 y];
		end
		zoom_level = zoom_level + 1;
		hist_time_range(zoom_level,:) = time_range;
		hist_dist_range(zoom_level,:) = dist_range;
	end
	if bot == 'o' % reset
		if zoom_level > 1
			zoom_level = zoom_level - 1;
		end
		time_range = hist_time_range(zoom_level,:);
		dist_range = hist_dist_range(zoom_level,:);
		isfill = 0;
	end
	if bot == 'O' % reset
		zoom_level = 1;
		hist_time_range(zoom_level,:) = ori_time_range;
		hist_dist_range(zoom_level,:) = ori_dist_range;
		time_range = hist_time_range(zoom_level,:);
		dist_range = hist_dist_range(zoom_level,:);
		isfill = 0;
		is_reduce_v = 0;
		single_norm = 1;
		amp = 5;
		ref_v = 10;
		is_dist = 1;
	end
	if bot == 'b'
		is_bin = ~is_bin;
	end
end

outevent.winpara = event.winpara;
outevent.isgood = event.isgood;
