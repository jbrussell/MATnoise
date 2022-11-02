% Run bootstrap for alpha and beta
% Randomly resample the data
%
% github.com/jbrussell
% 2022-03
%
% dlnbeta_dx: beta variation in longitude
% dlnbeta_dy: beta variation in latitude

function [alpha, dlnbeta_dx, dlnbeta_dy, beta_tau, azi_maxamp, alpha_err, dlnbeta_dx_err, dlnbeta_dy_err, beta_tau_err, azi_maxamp_err]=run_fit_alpha_beta_bootstrap_randresamp_gammacorr(azi,gamma,amp,amp_err,Nbs)

%     % Balanced resampling of all grid cells
%     % Each grid cell is selected an equal number of times
%     Ndata = length(amp);
%     xx = [1:Ndata]';
%     XX = repmat(xx,Nbs-1,1);
%     ind = randperm(Ndata*Nbs-1);
%     YY = XX(ind);
%     indxs = reshape(YY,Ndata,Nbs-1);
%     datamat = reshape(repmat(amp,1,Nbs-1),Ndata,Nbs-1);
%     AMP = [amp(:), datamat(indxs)];
%     datamat = reshape(repmat(azi,1,Nbs-1),Ndata,Nbs-1);
%     AZI = [azi(:), datamat(indxs)];
%     datamat = reshape(repmat(amp_err,1,Nbs-1),Ndata,Nbs-1);
%     AMP_err = [amp_err(:), datamat(indxs)];


    % Balanced resampling of events
    % Each event is sampled an equal number of times
    Nevs = size(amp,3);
    xx = [1:Nevs]';
    XX = repmat(xx,Nbs,1);
    ind = randperm(Nevs*Nbs);
    YY = XX(ind);
    indxs = reshape(YY,Nevs,Nbs);
    indxs = [xx, indxs];

    alpha_bs=zeros(Nbs,1);
    dlnbeta_dx_bs=zeros(Nbs,1);
    dlnbeta_dy_bs=zeros(Nbs,1);
    for ibs = 1:Nbs
        azi_bs = azi(:,:,indxs(:,ibs));
        gamma_bs = gamma(:,:,indxs(:,ibs));
        amp_bs = amp(:,:,indxs(:,ibs));
        amp_err_bs = amp_err(:,:,indxs(:,ibs));
        % Run fit on perturbed data
%         [~, alpha_bs(ibs), dlnbeta_dx_bs(ibs), dlnbeta_dy_bs(ibs)]=fit_alpha_beta(azi_bs(:),amp_bs(:),amp_err_bs(:));
        [~, alpha_bs(ibs), dlnbeta_dx_bs(ibs), dlnbeta_dy_bs(ibs)]=fit_alpha_beta_gammacorr(azi_bs(:),gamma_bs(:),amp_bs(:),amp_err_bs(:));
        
        if 0
            figure(1000); hold on;
            if ibs == 1
                clf;
            end
            plot(azi(:),amp(:),'xk');
            plot(azi_bs(:),amp_bs(:),'.')
            pause;
        end
    end
    beta_tau_bs = sqrt(dlnbeta_dx_bs.^2+dlnbeta_dy_bs.^2);
    azi_maxamp_bs = atan2d(dlnbeta_dx_bs,dlnbeta_dy_bs);
    
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
    
%     parastd=confint(para,.95);
%     dlnbeta_dx_err = (parastd(2,1)-parastd(1,1))/2;
%     dlnbeta_dy_err = (parastd(2,2)-parastd(1,2))/2;
%     beta_tau_err = ( (dlnbeta_dx.*dlnbeta_dx_err).^2 + (dlnbeta_dy.*dlnbeta_dy_err).^2 ).^0.5 ./ (dlnbeta_dx.^2+dlnbeta_dy.^2).^0.5;
%     azi_maxamp_err = 180/pi * ( (dlnbeta_dy.*dlnbeta_dx_err).^2 + (dlnbeta_dx.*dlnbeta_dy_err).^2 ).^0.5 ./ ( dlnbeta_dx.^2+dlnbeta_dy.^2 );
%     alpha_err = (parastd(2,3)-parastd(1,3))/2;

end

