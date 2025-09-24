function [tg, Ug, Ag, idxU, meta] = pick_ftan_ridges_robust(A, t, Tvec, dist, vmin, vmax, opts)
% PICK_FTAN_RIDGES_ROBUST
% Multi-branch FTAN ridge picker that operates in velocity domain.
%
% Inputs
%   A      : FTAN envelope map [Ntime x Nper] (non-negative)
%   t      : time vector (s) [Ntime x 1]
%   Tvec   : periods (s) [Nper x 1]
%   dist   : source?receiver distance (km)
%   vmin, vmax : group velocity bounds (km/s)
%   opts (struct, optional):
%       .Nu            (default 400)   % # velocity samples
%       .nBranches     (default 2)     % # dispersion branches to pick
%       .per_bounds       [nBranches x 2] vector defining period range allowed for each mode
%       .lambda_dU     (default 150)   % smoothness penalty on ?U (km/s)^2
%       .max_jump_U    (default 0.5)   % max |?U| per period step (km/s)
%       .w_amp         (default 1.0)   % weight on amplitude score
%       .normalizePerPeriod (true)     % robust z-score per period
%       .stopIfWeakZ   (default -Inf)  % stop if median z along path < thresh
%       .suppress_sigmaU (0.12)        % Gaussian half-width (km/s) to suppress around picked ridge
%       .suppress_gain (5)             % amount to subtract from Z along picked ridge
%
% Outputs
%   tg   : [Nper x K] arrival times (s) (tg = dist ./ Ug)
%   Ug   : [Nper x K] picked group velocities (km/s) for K branches
%   Ag   : [Nper x K] amplitudes of picks
%   idxU : [Nper x K] indices into U grid for each pick
%   meta : struct with fields Uvec, Z (z-scored map in velocity domain),
%          score arrays, and options used.
%
% Notes
%   - Builds a velocity-domain map by sampling A at t = dist ./ U.
%   - Uses DP with first-order smoothness and local jump constraint.
%   - Iteratively suppresses a narrow band around the picked ridge to find
%     additional branches (fundamental + overtones).
%
%
% jbrussell + ChatGPT - 9/2025

% ---------- Defaults ----------
if nargin < 7, opts = struct(); end
opts = set_defaults(opts, vmin, vmax, Tvec);

[Ntime, Nper] = size(A);
Uvec = linspace(vmin, vmax, opts.Nu).';           % [Nu x 1]
tau  = dist ./ Uvec;                              % [Nu x 1] times for each U

% ---------- Build velocity-domain amplitude map ----------
Avec = nan(opts.Nu, Nper);
for k = 1:Nper
    col = A(:,k);
    Avec(:,k) = interp1(t, col, tau, 'linear', NaN);
end

% ---------- Robust per-period z-score (stabilizes scoring) ----------
Z = Avec;
if opts.normalizePerPeriod
    for k = 1:Nper
        x = Avec(:,k);
        m = median(x, 'omitnan');
        s = median(abs(x - m), 'omitnan'); s = s + eps;
        Z(:,k) = (x - m)./s;
    end
else
    Z = log(Avec + eps);
end

% Invalid samples -> strong negative score so they don't get picked
Z(~isfinite(Z)) = -Inf;

% ---------- Multi-branch DP picking with suppression ----------
K = opts.nBranches;
Ug   = nan(Nper, K);
tg   = nan(Nper, K);
Ag   = nan(Nper, K);
idxU = nan(Nper, K);
path_scores = -inf(1,K);

Zwork = Z;

for b = 1:K
    [idx_path, best_score] = dp_path_velocity(Zwork, Uvec, opts);
    path_scores(b) = best_score;

    % Convert to outputs
    idxU(:,b) = idx_path;
    for k = 1:Nper
        ii = idx_path(k);
        % Skip if outside allowed period bounds
        if Tvec(k)<opts.per_bounds(b,1) || Tvec(k)>opts.per_bounds(b,2)
            ii = nan;
        end
        if ~isnan(ii)
            Ug(k,b) = Uvec(ii);
            tg(k,b) = dist / Ug(k,b);
            Ag(k,b) = Avec(ii,k);
        end
    end

    % Check strength
    z_along = arrayfun(@(k) Zwork(idxU(k,b),k), 1:Nper);
    if median(z_along(isfinite(z_along))) < opts.stopIfWeakZ
        % Too weak ? discard and stop
        Ug(:,b) = nan; tg(:,b) = nan; idxU(:,b) = nan; path_scores(b) = -Inf;
        break
    end

    % Suppress a Gaussian tube around the picked ridge in U
    Zwork = suppress_in_velocity(Zwork, Uvec, idxU(:,b), opts.suppress_sigmaU, opts.suppress_gain);
end

% ---------- Meta ----------
meta = struct('Uvec', Uvec, 'Z', Z, 'scores', path_scores, 'options', opts);

end

% ======================================================================
%                           Helper functions
% ======================================================================
function opts = set_defaults(opts, vmin, vmax, Tvec)
if ~isfield(opts,'Nu'),                opts.Nu = 400; end
if ~isfield(opts,'nBranches'),         opts.nBranches = 2; end
if ~isfield(opts,'per_bounds'),        opts.per_bounds = repmat([min(Tvec) max(Tvec)],2,1); end
if ~isfield(opts,'lambda_dU'),         opts.lambda_dU = 150; end
if ~isfield(opts,'max_jump_U'),        opts.max_jump_U = 0.5; end
if ~isfield(opts,'w_amp'),             opts.w_amp = 1.0; end
if ~isfield(opts,'normalizePerPeriod'),opts.normalizePerPeriod = true; end
if ~isfield(opts,'stopIfWeakZ'),       opts.stopIfWeakZ = -Inf; end
if ~isfield(opts,'suppress_sigmaU'),   opts.suppress_sigmaU = 0.12; end
if ~isfield(opts,'suppress_gain'),     opts.suppress_gain = 5; end
% Sanity clamp
opts.max_jump_U = max(0.05, min(opts.max_jump_U, max(0.8, (vmax-vmin)/2)));
end

function [idx_path, best_score] = dp_path_velocity(Z, Uvec, opts)
% Dynamic programming across velocity grid Uvec for each period.
[Nu, Nper] = size(Z);
dU = mean(diff(Uvec));
Jmax = max(1, ceil(opts.max_jump_U / dU)); % neighbors per step

dp = -inf(Nu, Nper);
bp = nan(Nu, Nper);

% Initialize
dp(:,1) = opts.w_amp * Z(:,1);

% Transitions
for k = 2:Nper
    for i = 1:Nu
        % restrict previous nodes to local neighborhood by ?U
        j1 = max(1, i - Jmax);
        j2 = min(Nu, i + Jmax);
        if j1 > j2, continue; end

        % cost from neighborhood
        dUi = (Uvec(i) - Uvec(j1:j2)).^2;
        cand_cost = dp(j1:j2, k-1) + opts.w_amp*Z(i,k) - opts.lambda_dU * dUi;

        [best, jrel] = max(cand_cost);
        dp(i,k) = best;
        bp(i,k) = j1 + jrel - 1;
    end
end

% Backtrack
[best_score, iend] = max(dp(:,end));
idx_path = nan(Nper,1);
idx_path(Nper) = iend;

for k = Nper:-1:2
    j = bp(idx_path(k), k);
    if isnan(j)
        [~, j] = max(dp(:,k-1));
    end
    idx_path(k-1) = j;
end
end

function Z2 = suppress_in_velocity(Z, Uvec, idx_path, sigmaU, gain)
% Subtract a Gaussian tube around the picked ridge in velocity space
Z2 = Z;
if ~isfinite(gain) || gain <= 0 || ~isfinite(sigmaU) || sigmaU <= 0
    return
end
[Nu, Nper] = size(Z);
for k = 1:Nper
    ii = idx_path(k);
    if isnan(ii), continue; end
    U0 = Uvec(ii);
    g = gain * exp(-0.5*((Uvec - U0)./sigmaU).^2);
    Z2(:,k) = Z2(:,k) - g;
end
end