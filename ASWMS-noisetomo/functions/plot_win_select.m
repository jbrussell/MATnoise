function plot_win_select(event,periods,winpara)
% This function is used to automatically select the window range used for gsdf method.
% The output format is
% v1 = winpara(1); t1 = winpara(2);
% v2 = winpara(3); t2 = winpara(4);
% and the window is defined by L/v1+t1 -- L/v2+t2

if ~exist('periods','var')
	setup_parameters;
	periods = parameters.periods;
end
if ~exist('winpara','var')
	winpara = event.winpara;
end

freqs = 1./periods;
minf = min(freqs);
maxf = max(freqs);

isdebug =1;

figure(1)
clf
hold on

amp = 0.4;
isgood = [event.stadata(:).isgood];
goodind = find(isgood > 0);
dist = [event.stadata(goodind).dist];
amp = (max(dist)-min(dist))/10*amp;
yrange = [ min(dist)-2*amp max(dist)+2*amp ];
trange = [ mean(dist)/5*0.8  mean(dist)/2*1.2 ];
trange = [ 0  mean(dist)/2*1.2 ];

for ista = 1:length(event.stadata)
    %  for ista = 20
    % set up time axis
	if event.stadata(ista).isgood < 0
		continue;
	end
    bgtime = event.stadata(ista).otime - event.otime;
    dt = event.stadata(ista).delta;
    Nt = length(event.stadata(ista).data);
    taxis = bgtime + [0:Nt-1]'*dt;
	data = event.stadata(ista).data;
    fN = 1/2/dt;
    [b,a] = butter(2,[minf/fN, maxf/fN]);
    data = filtfilt(b,a,data);
	data =  data./max(abs(data));
	plot(taxis,data*amp+event.stadata(ista).dist);
    
    [gausf,faxis] = build_gaus_filter(freqs,dt,Nt,0.06,0.1);
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
        norm_nband = abs(nband)./max(abs(nband));
        plot(taxis,norm_nband*amp+event.stadata(ista).dist,'k');
    end % end of loop ip
end % end of loop sta
if length(winpara)==4
plot([yrange/winpara(1)+winpara(2)],yrange,'r');
plot([yrange/winpara(3)+winpara(4)],yrange,'r');
end
plot([yrange/5],yrange,'k');
plot([yrange/2],yrange,'k');
ylim(yrange)
xlim(trange)
title(['Dist: ',num2str(mean(dist))]);
drawnow;
