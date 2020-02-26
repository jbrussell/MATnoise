function plot_SNR( ccf, grv_min, grv_max, r, ax )
% Calc SNR within the specified group velocity window
% JBR 7/8/18

snrdata = real(ifft((ccf)));
snrdata = fftshift(snrdata);
NN= length(snrdata);
lag = [-floor(NN/2):floor(NN/2)];
win_min = r./grv_max;
win_max = r./grv_min;
if win_min < 15; win_min = 0; end
if win_max < 50; win_max = 50; end
signal_ind = ((lag>-win_max & lag<-win_min) | (lag>win_min & lag<win_max));
signal_amp = sum(snrdata(signal_ind).^2)/length(snrdata(signal_ind));
noise_amp = sum(snrdata(~signal_ind).^2)/length(snrdata(~signal_ind));
snr = signal_amp/noise_amp;

plot(ax,lag,snrdata,'-','color',[0.5 0.5 0.5],'linewidth',2); hold on;
plot(ax,[win_min win_min],[-1 1]*max(abs(snrdata)),'-b','linewidth',2);
plot(ax,-[win_min win_min],[-1 1]*max(abs(snrdata)),'-b','linewidth',2);
plot(ax,[win_max win_max],[-1 1]*max(abs(snrdata)),'-b','linewidth',2);
plot(ax,-[win_max win_max],[-1 1]*max(abs(snrdata)),'-b','linewidth',2);
plot(ax,lag(signal_ind),snrdata(signal_ind),'-k','linewidth',2);
xlabel(ax,'Lag (s)');
title(ax,sprintf('SNR %.2f',snr));
xlim(ax,[-win_max win_max]*1.5);
set(ax,'linewidth',1.5,'fontsize',16);

end

