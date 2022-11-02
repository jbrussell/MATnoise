function [winpara outevent] = auto_win_select(event,mingroupv,maxgroupv,bandnum,center_freq)
% This function is used to automatically select the window range used for gsdf method.
% The output format is
% v1 = winpara(1); t1 = winpara(2);
% v2 = winpara(3); t2 = winpara(4);
% and the window is defined by L/v1+t1 -- L/v2+t2

isdebug = 0;

setup_parameters
setup_ErrorCode

periods = parameters.periods;
largest_epidist_range = parameters.largest_epidist_range;
is_removebadwindows = parameters.is_removebadwindows;

if ~exist('mingroupv')
    mingroupv = parameters.min_groupv;
end
if ~exist('maxgroupv')
    maxgroupv = parameters.max_groupv;
end
if ~exist('bandnum')
    bandnum = 20;
end
if ~exist('center_freq')
    center_freq = parameters.cent_freq;
end

cycle_before = parameters.cycle_before;
cycle_after = parameters.cycle_after;
min_dist_tol = parameters.min_dist_tol;
max_dist_tol = parameters.max_dist_tol;

peakamptol = 0.5;   % the ratio of accepted peak compare to the largest peak
peak_search_range = 1;
groupv_diff_tol = 0.2;
positive_disp_weight = 5;
min_sta_num = parameters.min_sta_num; %10;
bad_f_num_tol = bandnum*0.6;

minf = 1/periods(end);
maxf = 1/periods(1);
freqs = linspace(minf,maxf,bandnum);
[temp center_freq_index] = min(abs(freqs - center_freq));


% set up some useful arrays
stlas = [event.stadata(:).stla];
stlos = [event.stadata(:).stlo];
stnms = {event.stadata(:).stnm};
dists = [event.stadata(:).dist];

% check whether the stations are in the range
for ista = 1:length(event.stadata)
	event.stadata(ista).isgood = 1;
	if ~Is_inrange(stlas(ista),stlos(ista),parameters)
		event.stadata(ista).isgood = ErrorCode.sta_outofrange;
	end
end
    
% if the array is too large and sparse, select the major part of it.
isgood = [event.stadata(:).isgood];
goodind = find(isgood>0);
dist = [event.stadata(goodind).dist];
if mean(dist) < min_dist_tol || mean(dist) > max_dist_tol
    disp(['Event: ',event.id,' is not in the proporal range']);
    winpara = 0;
    outevent = event;
    return;
end
if max(dist) - min(dist) > largest_epidist_range
	disp(['Too large array size, only use part of stations for window picking']);
	distbin = min(dist):largest_epidist_range/10:max(dist);
	[stanum,middist] = hist(dist,distbin);
	% find the epicenter distance range contains most stations.
	maxstanum = 0;
	for ibin = 1:length(middist) - 10
		sumstanum = sum(stanum(ibin:ibin+10));
		if sumstanum > maxstanum
			maxstanum = sumstanum;
			maxibin = ibin;
		end
	end
	distrange = [middist(maxibin)-largest_epidist_range/20,middist(maxibin+10)+largest_epidist_range/20];
	dist = [event.stadata(:).dist];
	outind = find(dist < distrange(1) | dist > distrange(2));
	for iout = outind
		event.stadata(iout).isgood = ErrorCode.sta_outofepidist;
	end
end

good_sta_num = 0;
for ista = 1:length(event.stadata)
    %  for ista = 20
    % set up time axis
    bgtime = event.stadata(ista).otime - event.otime;
    dt = event.stadata(ista).delta;
    Nt = length(event.stadata(ista).data);
    taxis = bgtime + [0:Nt-1]'*dt;
    
    % Build up gaussian filter
    [gausf,faxis] = build_gaus_filter(freqs,dt,Nt,parameters.min_width,parameters.max_width);
    
    % get original data and make the fourier transform
    odata = event.stadata(ista).data;
    if size(odata,1) == 1  % in matlab, the fast direction is column
        odata = odata';
    end
    fftodata = fft(odata);
    
    clear envelop_nbands nband nbands norm_envelop
    % apply narrow-band filters
    for ip = 1:length(freqs);
        nband = fftodata .* [gausf(:,ip); zeros(Nt-length(gausf(:,ip)),1)];
        nband = ifft(nband);
        nbands(:,ip) = nband;
    end % end of loop ip
    
    % remove the area out of window defined by min and max group velocity
    tmin = event.stadata(ista).dist/maxgroupv;
    tmax = event.stadata(ista).dist/mingroupv;
    if tmax > 5500
        tmax = 5500;
    end
    if tmin < taxis(1) || tmax > taxis(end)
        disp(['Station ',event.stadata(ista).stnm,' does not contain enough data']);
        if ista == length(event.stadata)
            groupdelay(:,ista) = 0;
            snr(:,ista) = 0;
        end
        continue;
    end
	if event.stadata(ista).isgood < 0
        if ista == length(event.stadata)
            groupdelay(:,ista) = 0;
            snr(:,ista) = 0;
        end
        continue;
	end
    good_sta_num = good_sta_num+1;
    inwin_ind = find(taxis > tmin & taxis < tmax);
    outwin_ind = find(taxis < tmin | taxis > tmax);
    envelop_nbands=abs(nbands);
    envelop_nbands(outwin_ind,:) = 0; % apply group velocity window
    for ip = 1:length(freqs) % Normalize it for grading purpose
        envelop_nbands(:,ip) = envelop_nbands(:,ip) / max(envelop_nbands(:,ip));
    end
    
    
    % Find local maximum as potential peaks
    clear peaks
    for ip = 1:length(freqs)
        envelopfun = envelop_nbands(:,ip);
        diff_fun = diff(envelopfun);
        peakind = inwin_ind(find(diff_fun(inwin_ind-1)>0 & diff_fun(inwin_ind)<0));
        bigpeakamp = max(envelopfun);
        peakind = peakind(find(envelopfun(peakind)>bigpeakamp*peakamptol));
        peaks(ip).peaktimes = taxis(peakind);
        peaks(ip).peakamps = envelopfun(peakind);
        peaks(ip).peaknum = length(peakind);
    end
    
    % From center frequency band, search for the best peak to both sides
    ip = center_freq_index;
    [temp ind] = max(peaks(ip).peakamps);
    peaks(ip).bestpeak = peaks(ip).peaktimes(ind);
    peaks(ip).bestpeaksnr = peaks(ip).peakamps(ind)/sum([peaks(ip).peakamps]);
    % To the lower frequency bands
    for ip=center_freq_index-1:-1:1  % loop for lower frequencies
        if isempty(peaks(ip).peaktimes) ...
            || isempty(peaks(ip+1).peaktimes) || isempty(peaks(ip+1).bestpeak)
            peaks(ip).bestpeak = peaks(ip+1).bestpeak;
            peaks(ip).bestpeaksnr = 0;
            continue;
        end
        [temp closest_peak_i] = min(abs(peaks(ip+1).bestpeak - peaks(ip).peaktimes));
        bestpeaki = closest_peak_i;
        bestpeakamp = peaks(ip).peakamps(bestpeaki);
        groupV0 = event.stadata(ista).dist / peaks(ip+1).bestpeak;
        for ipeak = closest_peak_i-peak_search_range:closest_peak_i+peak_search_range
            if ipeak >= 1 && ipeak <= length(peaks(ip).peaktimes) % search valid
                groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(ipeak);
                if groupvi > groupV0*(1-groupv_diff_tol) && ...
                        groupvi < groupV0*(1+groupv_diff_tol)
                    if peaks(ip).peakamps(ipeak) > bestpeakamp
                        bestpeaki = ipeak;
                        bestpeakamp = peaks(ip).peakamps(ipeak);
                    end % end for finding largest nearby peak
                end % end for check group v
            end % end for check ipeak vavid
        end % end of loop ipeak
        peaks(ip).bestpeak = peaks(ip).peaktimes(bestpeaki);
        peaks(ip).bestpeaksnr = peaks(ip).peakamps(bestpeaki)/sum([peaks(ip).peakamps]);
        groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(bestpeaki);
        if groupvi < groupV0*(1-groupv_diff_tol) || ...
                groupvi > groupV0*(1+groupv_diff_tol)
            peaks(ip).bestpeak = peaks(ip+1).bestpeak;
            peaks(ip).bestpeaksnr = 0;
        end
    end % end of loop ip
    % To the higher frequency bands
    for ip=center_freq_index+1:length(freqs)  % loop for higher frequencies
        if isempty(peaks(ip).peaktimes) ...
                || isempty(peaks(ip-1).peaktimes) || isempty(peaks(ip-1).bestpeak)
            peaks(ip).bestpeak = peaks(ip-1).bestpeak;
            peaks(ip).bestpeaksnr = 0;
            continue;
        end
%         disp([num2str(ista),',',num2str(ip)]);
        [temp closest_peak_i] = min(abs(peaks(ip-1).bestpeak - peaks(ip).peaktimes));
        bestpeaki = closest_peak_i;
        amp_point = peaks(ip).peakamps(bestpeaki);
        groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(bestpeaki);
        groupV0 = event.stadata(ista).dist / peaks(ip-1).bestpeak;
        if groupvi/groupV0 > 1
            v_point = 1-(groupvi/groupV0-1)/groupv_diff_tol;
        else
            % double the weight for having a positive dispersion curve
            v_point = 1-(1-groupvi/groupV0)/groupv_diff_tol/positive_disp_weight;
        end
        best_sum_point = amp_point + v_point;
        for ipeak = closest_peak_i-peak_search_range:closest_peak_i+peak_search_range
            if ipeak >= 1 && ipeak <= length(peaks(ip).peaktimes) % search valid
                groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(ipeak);
                if groupvi > groupV0*(1-groupv_diff_tol) && ...
                        groupvi < groupV0*(1+groupv_diff_tol)
                    amp_point = peaks(ip).peakamps(ipeak);
                    if groupvi/groupV0 > 1
                        v_point = 1-(groupvi/groupV0-1)*10;
                    else
                        % double the weight for having a positive dispersion curve
                        v_point = 1-(1-groupvi/groupV0)*10/positive_disp_weight;
                    end
                    %                     disp([num2str(ip),':',num2str(amp_point),',',num2str(v_point)]);
                    if amp_point + v_point > best_sum_point
                        bestpeaki = ipeak;
                        best_sum_point = amp_point + v_point;
                    end % end for finding largest nearby peak
                end % end for check group v
            end % end for check ipeak vavid
        end % end of loop ipeak
        peaks(ip).bestpeak = peaks(ip).peaktimes(bestpeaki);
        peaks(ip).bestpeaksnr = peaks(ip).peakamps(bestpeaki)/sum([peaks(ip).peakamps]);
        groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(bestpeaki);
        if groupvi < groupV0*(1-groupv_diff_tol) || ...
                groupvi > groupV0*(1+groupv_diff_tol)
            peaks(ip).bestpeak = peaks(ip-1).bestpeak;
            peaks(ip).bestpeaksnr = 0;
        end
    end % end of loop ip
    
    for ip=1:length(freqs)
        if isempty(peaks(ip).bestpeak)
            peaks(ip).bestpeak = 0;
            peaks(ip).bestpeaksnr = 0;
        end
    end

    % record the group delay
    groupdelay(:,ista) = [peaks(:).bestpeak]';
    snr(:,ista) = [peaks(:).bestpeaksnr]';
    if isdebug
        stadata(ista).peaks = peaks;
        stadata(ista).envelop_nbands = envelop_nbands;
    end
end % end of loop sta

if good_sta_num < min_sta_num
    disp(['Event: ',event.id, ' doesn''t have enough stations, skip!']);
    winpara =0;
	outevent = event;
    return
end

dist = [event.stadata(:).dist];
for ip = 1:length(freqs)
    [groupv(ip) offset(ip)] = groupv_fit(dist,groupdelay(ip,:),snr(ip,:),mingroupv,maxgroupv);
end % end of ip

clear bgtime endtime
for ista = 1:length(event.stadata)
    peaktimes = dist(ista)./groupv + offset;
    bgtimes = peaktimes - cycle_before./freqs;
    endtimes = peaktimes + cycle_after./freqs;
    bgtime(ista) = min(bgtimes);
    endtime(ista) = max(endtimes);
end

if isdebug
    figure(39)
    clf
    hold on
    plot(dist,bgtime,'o');
    plot(dist,endtime,'ro');
end
isgood = [event.stadata(:).isgood];
goodind = find(isgood > 0);
para = polyfit(dist(goodind),bgtime(goodind),1);
v1 = 1/para(1);
t1 = para(2);
para = polyfit(dist(goodind),endtime(goodind),1);
v2 = 1/para(1);
t2 = para(2);
winpara = [v1,t1,v2,t2];
if 0
    for ista = 1:length(event.stadata)
        envelop_nbands = stadata(ista).envelop_nbands;
        peaks = stadata(ista).peaks;
        if isempty(envelop_nbands)
            continue;
        end
        for ip = 1:length(freqs)
            norm_envelop(:,ip) = envelop_nbands(:,ip) / max(envelop_nbands(:,ip));
        end
        figure(36)
        clf
        hold on
        [xi yi] = ndgrid(taxis,freqs);
        contourf(xi,yi,norm_envelop);
        for ip = 1:length(freqs)
            plot(peaks(ip).bestpeak,freqs(ip),'gx','markersize',15,'linewidth',3);
        end
        for ip = 1:length(freqs)
            plot(peaks(ip).peaktimes,ones(size(peaks(ip).peaktimes))*freqs(ip),'r.','markersize',15);
        end
        bgt = dist(ista)/v1+t1;
        endt = dist(ista)/v2+t2;
        plot([bgt bgt],[freqs(1) freqs(end)],'r','linewidth',3);
        plot([endt endt],[freqs(1) freqs(end)],'r','linewidth',3);
        
        xlim([dist(ista)/maxgroupv, dist(ista)/mingroupv])
        pause
    end
end
if 0
    for ip=1:length(freqs)
        figure(38)
        clf
        hold on
        plot(dist,groupdelay(ip,:),'x');
        plot([min(dist),max(dist)],...
            [min(dist)/groupv(ip)+offset(ip),max(dist)/groupv(ip)+offset(ip)])
        title(['V:',num2str(groupv(ip)),' t:',num2str(offset(ip))]);
        pause
    end
end

if is_removebadwindows==1
    bad_f_ind = find(groupv == mingroupv | groupv == maxgroupv);
    if length(bad_f_ind) > bad_f_num_tol
        disp(['Event: ',event.id, ' SNR is too low, skip!']);
        winpara = 0;
        outevent = event;
        return
    end
end

outevent = event;

end % end of function
