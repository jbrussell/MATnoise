function [ snr ] = calc_SNR( sac, parameters, isfigure_snr )
% Calc SNR within the specified group velocity window
% JBR 9/28/19
prefilter = parameters.prefilter;
refv = parameters.refv;
data = sac.DATA1;
dt = sac.DELTA;
NN= length(data);
t = sac.B + [0:NN-1]*sac.DELTA;

data = data.*tukeywin(NN,0.08);
fN = 1/2/dt;
[b,a] = butter(2,[1/prefilter(2)/fN, 1/prefilter(1)/fN]);
data = filtfilt(b,a,data);

grv_min = refv*0.5;
grv_max = refv*1.1;
r = vdist(sac.STLA,sac.STLO,sac.EVLA,sac.EVLO)/1e3;
win_min = r./grv_max;
win_max = r./grv_min;
signal_ind = (t>win_min & t<win_max);
signal_amp = sum(data(signal_ind).^2)/length(data(signal_ind));
noise_amp = sum(data(~signal_ind).^2)/length(data(~signal_ind));
snr = signal_amp/noise_amp;

if isfigure_snr
    figure(101); clf;
    subplot(2,1,1);
    plot(t,data,'-','color',[0.5 0.5 0.5],'linewidth',1); hold on;
    plot([win_min win_min],[-1 1]*max(abs(data)),'-r','linewidth',2);
    plot([win_max win_max],[-1 1]*max(abs(data)),'-r','linewidth',2);
    plot(t(signal_ind),data(signal_ind),'-k','linewidth',1);
    xlabel('Lag (s)');
    title(sprintf('SNR %.2f',snr));
    xlim([win_min-500 win_max+500]);
    set(gca,'linewidth',1.5,'fontsize',16);

    subplot(2,1,2);
    plot(t,data,'-','color',[0.5 0.5 0.5],'linewidth',1); hold on;
    plot([win_min win_min],[-1 1]*max(abs(data)),'-r','linewidth',2);
    plot([win_max win_max],[-1 1]*max(abs(data)),'-r','linewidth',2);
    plot(t(signal_ind),data(signal_ind),'-k','linewidth',1);
    xlabel('Lag (s)');
    xlim([min(t) max(t)]);
    set(gca,'linewidth',1.5,'fontsize',16);
    
    pause;
end

end

