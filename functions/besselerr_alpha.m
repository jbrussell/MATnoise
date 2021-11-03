function err = besselerr_alpha(alph,c,xsp,omega,r)
% Error function for fitting average bessel function with distance for a 
% single frequency in order to determine attenuation parameter, alpha. This
% should be run after phase velocity has been estimated at each frequency!
%
% Re(rho) = J_0(omega*r/c) * exp(-alpha*r)
%
% Modified on 10/2021 by JBR
%

Isfigure=0;

% x1 = omega.*tw;
x1 = omega .* r ./ c;

% First part of error: fit the Bessel Function
%F1 = normalise(xsp);
F1 = xsp;
F_env = abs(hilbert(F1));
% F_env = SmoothAnalyticEnv(r,F1);
log_F_env = log(F_env);

be = besselj(0,x1) .* exp(-alph.*r);
be = be./mean(abs(be)).*mean([abs(F1)]);
be_env = abs(hilbert(be));
% be_env = SmoothAnalyticEnv(r,be);
log_be_env = log(be_env);
% if is_normbessel
% %     be = be./abs(hilbert(be));
%     be = be./SmoothAnalyticEnv(waxis/2/pi,be);
% end
F1z =  log_be_env(:) - log_F_env(:);
% F1z =  be_env(:) - F_env(:);


% % ORIGINAL
err = F1z(:);

if 1
    figure(99); clf;
    subplot(2,1,1);
    hold on;
    plot(r,F1,'.k');
    plot(r,be,'.y');
    plot(r,F_env,'.r');
    plot(r,SmoothAnalyticEnv(r,F1),'.m');
    plot(r,be_env,'.b');
    plot(r,SmoothAnalyticEnv(r,be),'.c');
    
    subplot(2,1,2); hold on;
    plot(r,log_F_env,'.r');
    plot(r,log_be_env,'.b');
end
    

if Isfigure>0
    figure(Isfigure)
    clf
    subplot(2,1,1)
    hold on
    plot(r,log_F_env,'.','linewidth',2);
    plot(r,log_be_env,'.','linewidth',2);
    
    pause;
end
end
