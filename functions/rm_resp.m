function outdata = rm_resp(data,zeros,poles,gain,delta)
%% function to remove instrument response given poles and zeros
%
% INPUT:
% data: raw uncorrected data
% zeros: zeros of transfer function
% poles: poles of transfer function
% gain: Gain defined by sensitivity * A0
% delta: 1/sampleRate [s]
%
% OUTPUT:
% outdata: Data after removal of instrument response
%
% written by Ge Jin, 2014/02/27
% Modified by jbrussell: 2020/10/30
%
% intrace = intrace(1);
isfigure = 0;
lo_corner = 0.005;  % in Hz
npoles=5; 

% data = intrace.data;
data = detrend(data);
data = flat_hanning_win(1:length(data),data,1,length(data),50);

% T = (intrace.endTime - intrace.startTime)*24*3600;
N = length(data);
T = N*delta; %[0:N-1]*delta
% delta = 1/intrace.sampleRate;

if mod(N,2)
     faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
else
     faxis = [0:N/2,-N/2+1:-1]*(1/T);
end

% poles = intrace.sacpz.poles;
% zeros = intrace.sacpz.zeros;
% gain = intrace.sacpz.constant;
w = faxis.*2*pi;
resp = ones(size(w));
for ip = 1:length(poles)
	resp = resp./(i*w - poles(ip));
end
for ip = 1:length(zeros)
	resp = resp.*(i*w - zeros(ip));
end
resp = resp*gain;
if isfigure
	figure(33)
	clf
	set(gcf,'position',[360   514   900   400]);
	hold on
	subplot(1,2,1)
	set(gca,'fontsize',18)
	semilogy(faxis,abs(resp),'rx');
	subplot(1,2,2)
	set(gca,'fontsize',18)
	plot(faxis,angle(resp),'rx');
end

lo_w=2*pi*lo_corner;
hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
norm_trans=hpfiltfrq./resp;    % this is normalization transfer function
norm_trans(find(isnan(norm_trans))) = 0;

fftdata = fft(data);
fftdata = fftdata(:).*norm_trans(:);
data_cor = real(ifft(fftdata));

outdata = data_cor;

% disp(['Station: ',intrace.station,'.',intrace.channel,' deconv to ',intrace.sacpz.units]);

return

