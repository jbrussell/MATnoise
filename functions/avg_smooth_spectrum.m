function [ power_avg_sm ] = avg_smooth_spectrum( R, Z, dt, sm_type )
% Calculate the average smooth spectrum given R and Z data
%
% Josh Russell
% github.com/jbrussell

N = length(R);
f = [0:(N-mod(N-1,2))/2 , -(N-mod(N,2))/2:-1]/dt/N;

fftR = fft(R);
fftZ = fft(Z);

% Index positive values (single-sided)
I = find(f>=0);
f = f(I);
fftR = fftR(I);
fftZ = fftZ(I);

if strcmp(sm_type,'octavesmooth')
    Noct_sm = 10; % 1/N octave smoothing (smaller number means more smoothing)
    power_R_sm = smoothSpectrum(abs(fftR(:)),f(:),Noct_sm);
    power_Z_sm = smoothSpectrum(abs(fftZ(:)),f(:),Noct_sm);
elseif strcmp(sm_type,'defaultsmooth')
    npts = 10/dt; %N/100;
    power_R_sm = smooth(abs(fftR(:)),npts);
    power_Z_sm = smooth(abs(fftZ(:)),npts);
end

power_avg_sm = 0.5*(power_R_sm + power_Z_sm);
power_avg_sm = power_avg_sm(:);

if 0
    figure(1); clf;
    t = [0:length(R)-1]*dt;

    subplot(2,1,1); box on; hold on;
    plot(t,R,'-r');
    plot(t,Z,'-b');

    subplot(2,1,2);  box on; hold on;
    plot(f,abs(fftR),'-r');
    plot(f,abs(fftZ),'-b');
    plot(f,power_R_sm,'-m','linewidth',2);
    plot(f,power_Z_sm,'-c','linewidth',2);
    plot(f,power_avg_sm,'-k','linewidth',2);
    set(gca,'yscale','log','xscale','log');
end

% Convert back to double-sided (0,+,-)
power_avg_sm = [power_avg_sm(:); flip(power_avg_sm(2:end))]';

end

