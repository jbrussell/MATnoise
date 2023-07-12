function [wavefield] = microseism_source(t,r,vel,freq,phi_rand)
% Microseism source from equation (4) of S. Schippkus et al. (2023)
%
% t: time axis [s]
% r: distance from source [km]
% vel: velocity of medium [km/s]
% freq: 1xN vector of frequencies to sum over [Hz]
% phi_rand: 1xN vector of random uniform phase between [0,2*pi] rad
%
% jbrussell - 7/2023

Nt = length(t);
Nfreq = length(freq);

% Adjust sizes
freq = freq(:);
freq = repmat(freq,1,Nt);
phi_rand = phi_rand(:);
phi_rand = repmat(phi_rand,1,Nt);
tvec = t(:)';
t = repmat(t,Nfreq,1);

% Green's function (delta function)
gf = zeros(size(tvec));
t_arrival = r / vel;
[~,I] = min(abs(tvec-t_arrival));
gf(I) = 1;

% Calculate microseism source (Gualtieri et al.; 2020)
A = ones(size(freq)); % assume uniform
source = sum(A .* cos(2*pi*freq.*t + phi_rand),1);

% Convolve Green's function with source to get seismogram
% wavefield = conv(gf,source,'same');
wavefield = real(ifft(fft(gf).*fft(source)));

% Plot example
if 0
    figure(999); clf;
    set(gcf,'position',[616   199   649   819]);
    
    subplot(4,1,1); box on; hold on;
    plot(tvec,gf,'b','linewidth',2);
    plot(t_arrival,[0 1],'r','linewidth',1);
    set(gca,'fontsize',15,'linewidth',1.5);
    title('$G(x,x_N) = \delta(t-L/c)$','Interpreter','latex','fontsize',18)
    
    subplot(4,1,2); box on; hold on;
    plot(tvec,source,'b','linewidth',2);
    plot(t_arrival,[min(source) max(source)],'r','linewidth',1);
    set(gca,'fontsize',15,'linewidth',1.5);
    title('$S_N(t,x_N) = \sum_{i}\cos(2\pi f_i t + \phi_i$)','Interpreter','latex','fontsize',18)
    
    subplot(4,1,3); box on; hold on;
    plot(tvec,wavefield,'b','linewidth',2);
    plot(t_arrival,[min(wavefield) max(wavefield)],'r','linewidth',1);
    set(gca,'fontsize',15,'linewidth',1.5);
    title('$W(t,x,x_N) = G(t,x,x_N) * S_N(t,x_N)$','Interpreter','latex','fontsize',18)
    
    subplot(4,1,4); box on; hold on;
    plot(tvec,ifft(fft(wavefield).*conj(fft(source))),'b','linewidth',2);
    plot(t_arrival,[0 1],'r','linewidth',1);
    set(gca,'fontsize',15,'linewidth',1.5);
    xlabel('Time (sec)');
    title('$W(t,x,x_N) \star S_N(t,x_N)$','Interpreter','latex','fontsize',18)
    
%     pause
    drawnow
end

end

