function [ data_filt ] = tukey_filt( data,coperiod,dt,costap_wid )
% JBR -- 10/21/16
% Apply cosine taper (Tukey window) in the frequency domain
%
% costap_wid: 0 => square window
%             1 => Hann window (looks similar to Gaussian)
%             Default is 0.5
%


% BUILD TUKEYWIN FILTER
Nt = length(data);
faxis = [0:1/Nt:1/dt/2,-1/dt/2+1/Nt:1/Nt:-1/Nt];
fmax = 1/coperiod(1);
fmin = 1/coperiod(2);
[~ , If] = find((faxis>=fmin & faxis<=fmax) | (faxis<=-fmin & faxis>=-fmax));
costap_len = length(If);
costap_full = zeros(size(data));
costap = tukeywin(costap_len,costap_wid);
costap_full(If) = costap;

% APPLY FILTER
data_filt = data.*costap_full;

if 0
    figure(99); clf;
    plot(1./faxis(1:Nt-1),abs(data(1:Nt-1))','-k'); hold on;
    plot(1./faxis(1:Nt-1),abs(data_filt(1:Nt-1))','-r');
    xlabel('Period');
    pause;
end


end

