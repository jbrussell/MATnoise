% This function is used for fitting the azimuthal dependence of the
% corrected amplitude decay term to extract the attenuation coefficent
% (alpha) and the local site amplification (beta)
%
% Equation 9 in Bao et al. (2016) GJI; doi:10.1093/gji/ggw151
%
% dlnbeta_dx*sind(azi) + dlnbeta_dy*cosd(azi) - alpha = amp_term
%
% github.com/jbrussell
% 2022-03
%
% dlnbeta_dx: beta variation in longitude
% dlnbeta_dy: beta variation in latitude

function [alpha, dlnbeta_dx, dlnbeta_dy, dbeta_beta, azi_maxamp, alpha_std, dlnbeta_dx_std, dlnbeta_dy_std, dbeta_beta_std, azi_maxamp_std] = fit_alpha_beta_lsqr(azi,amp,amp_err)
    azi = azi(:);
    amp = amp(:);
    amp_err = amp_err(:);
    
    Nobs = length(azi);
    G = zeros(Nobs,3);
    for ii = 1:Nobs
        G(ii,:) = [sind(azi(ii)) cosd(azi(ii)) -1];
    end
    % remove nan values;
    inan = find(isnan(amp) | isnan(amp_err));
    % igood = find(~isnan(dlnbeta));
    G(inan,:) = [];
    amp(inan,:) = [];
    amp_err(inan,:) = [];
    azi(inan,:) = [];
    
    % Weighting by uncertainties
    W = diag(1./amp_err);
%     W = eye(length(amp_err));
    
    % small amount of damping
    dampweight0 = 0;
    F = eye(size(G,2));
    f = zeros(size(G,2),1);
    % Rescale the kernel
    NR=norm(F,1);
    NA=norm(W*G,1);
    dampweight = dampweight0*NA/NR;
    
    % Construct full G matrix
    H = [W*G; dampweight*F];
    h = [W*amp; dampweight*f];
%     H = [W*G];
%     h = [W*amp];
    
    % Invert for model parameters
    params = (H'*H)\H'*h;
    dlnbeta_dx = params(1);
    dlnbeta_dy = params(2);
    alpha = params(3);
    
    % Calculate model uncertainties
%     Ginv = (H'*H)\G'*W;
%     params_std = diag(Ginv*diag(amp_err.^2)*Ginv').^(1/2);
    params_std = diag(inv(H'*H)).^(1/2);
%     params_std = diag(Ginv*Ginv').^(1/2);
    dlnbeta_dx_std = params_std(1);
    dlnbeta_dy_std = params_std(2);
    alpha_std = params_std(3);

    % Lin et al. (2012) equivalent
    dbeta_beta = sqrt(dlnbeta_dx.^2+dlnbeta_dy.^2);
    azi_maxamp = atan2d(dlnbeta_dx,dlnbeta_dy);
    
    % Transform errors
    dbeta_beta_std = ( (dlnbeta_dx.*dlnbeta_dx_std).^2 + (dlnbeta_dy.*dlnbeta_dy_std).^2 ).^0.5 ./ (dlnbeta_dx.^2+dlnbeta_dy.^2).^0.5;
    azi_maxamp_std = 180/pi * ( (dlnbeta_dy.*dlnbeta_dx_std).^2 + (dlnbeta_dx.*dlnbeta_dy_std).^2 ).^0.5 ./ ( dlnbeta_dx.^2+dlnbeta_dy.^2 );

% % Compare uncertainty estimates with built in Matlab Function    
%     [x,stdx_scaled,mse,S_scaled] = lscov(G,amp,diag(amp_err.^2));
%     stdx = stdx_scaled * sqrt(1/mse);
	
%     % Estimate chi2 misfit
%     dlnbeta_pre = G * params;
%     e = (dlnbeta - dlnbeta_pre) ./ dlnbeta_err;
%     chi2 = (e'*e)/length(dlnbeta);
%     
%     % Calculate model resolution and chi2
%     Ginv = (A'*A)\mat';
%     R = Ginv * mat; % model resolution
%     D = mat * Ginv; % data resolution
%     % Effective degrees of freedom
%     v = length(dt) - trace(D);
% %         v = trace(D);
%     % normalized chi2 uncertainties
%     res = (dt-mat*phaseg);
%     res(diag(W)==0) = nan;
%     rms_res = sqrt(nanmean(res.^2));
%     dt_std = rms_res;
%     chi2 = nansum(res.^2./dt_std.^2)/v;


    if 0
        figure(100); clf;
        plot(azi,amp,'xb'); hold on;
        x = [0:360];
        plot(x,dlnbeta_dx*sind(x) + dlnbeta_dy*cosd(x) - alpha,'-r');
    end

end
