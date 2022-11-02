% Run bootstrap for alpha and beta 
% Perturb data by randomly sampling from a Gaussian distribution centered
% on the observed value with standard deviation equal to the uncertainty
%
% github.com/jbrussell
% 2022-03
%
% dlnbeta_dx: beta variation in longitude
% dlnbeta_dy: beta variation in latitude

function [alpha, dlnbeta_dx, dlnbeta_dy, beta_tau, azi_maxamp, alpha_err, dlnbeta_dx_err, dlnbeta_dy_err, beta_tau_err, azi_maxamp_err]=run_fit_alpha_beta_bootstrap(azi,amp,amp_err,Nbs)
    
    alpha_bs = [];
    dlnbeta_dx_bs = [];
    dlnbeta_dy_bs = [];
    beta_tau_bs = [];
    azi_maxamp_bs = [];
    for ibs = 1:Nbs
        %perturb data
        amp_bs = normrnd(amp, amp_err);
        
        [~, alpha_bs(ibs), dlnbeta_dx_bs(ibs), dlnbeta_dy_bs(ibs)]=fit_alpha_beta(azi(:),amp_bs(:),amp_err(:));
    end
    beta_tau_bs = sqrt(dlnbeta_dx_bs(ibs).^2+dlnbeta_dy_bs(ibs).^2);
    azi_maxamp_bs = atan2d(dlnbeta_dx_bs(ibs),dlnbeta_dy_bs(ibs));
    
    alpha = median(alpha_bs);
    dlnbeta_dx = median(dlnbeta_dx_bs);
    dlnbeta_dy = median(dlnbeta_dy_bs);
    beta_tau = median(beta_tau_bs);
    azi_maxamp = median(azi_maxamp_bs);
    alpha_err = std(alpha_bs);
    dlnbeta_dx_err = std(dlnbeta_dx_bs);
    dlnbeta_dy_err = std(dlnbeta_dy_bs);
    beta_tau_err = std(beta_tau_bs);
    azi_maxamp_err = std(azi_maxamp_bs);

end

