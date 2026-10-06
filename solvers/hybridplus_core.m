function result = hybridplus_core(problem, x0, options, variant)
%HYBRIDPLUS_CORE Shared implementation of restarted hybridplus.
%
%   RESULT = HYBRIDPLUS_CORE(PROBLEM, X0, OPTIONS, VARIANT) with VARIANT:
%     'restarted'  Restarted hybridplus (Algorithm 1 of the manuscript).
%     'identity'   hybridplus_cr: restarted hybridplus with B_k = I for all
%                  k (ablation of the diagonal scaling).
%
%   Method (manuscript Section 3). With s = x_k - x_{k-1}, y = F_k - F_{k-1},
%     b_i = max(ell, min(y_i/s_i, u)) if s_i ~= 0, b_i = 1 otherwise,
%     B_k = diag(1./b),   B_0 = I,   d_0 = -F_0,
%     p = d_{k-1},  w = y - (p'y/||p||^2) p,
%     z = p + sqrt(1-t^2)*||p||/(t*||w||) * w,
%     beta = max(0, F_k' B_k z) / max(p'z, ||F_{k-1}||^2),
%     d_k = -B_k F_k + beta*p.
%   Restart (beta = 0, d_k = -B_k F_k) iff p == 0, or w == 0, or the
%   computed d_k == 0. All three tests are EXACT zero tests, as in the
%   convergence theory.
%
%   Line search: from lambda = lambda0_k, test x_k + lambda d_k, then
%   x_k - lambda d_k; accept the first trial with
%     f(trial) <= (1 + tau_k) f(x_k) - sigma lambda^2 ||d_k||^2,
%   f = 0.5||F||^2; otherwise lambda <- rho*lambda. The memory uses the
%   unsigned direction d_k even when the negative trial is accepted.
%
%   Initial trial step lambda0_k (option initialStep; manuscript eq. for
%   the safeguarded initial trial). Every rule is clipped to
%   [stepMin, stepMax] with stepMin <= 1 <= stepMax; a zero, non-finite or
%   undefined estimate falls back to 1. The convergence theory holds for
%   any rule obeying these bounds.
%     'unit'    lambda0_k = 1.
%     'fd'      lambda0_k = |F_k'd_k| / |d_k'(F(x_k + eps d_k) - F_k)/eps|,
%               eps = probeEpsilon: a finite-difference estimate of
%               |F'd|/|d'J d| (the initial trial of Li and Li 2011).
%               Costs ONE counted residual evaluation per iteration.
%     'secant'  lambda0_k = |F_k'd_k| / (mu_k ||d_k||^2), with the secant
%               curvature mu_k = |s_{k-1}'y_{k-1}| / ||s_{k-1}||^2 along the
%               last step; lambda0_0 = 1. No extra evaluation.
%
%   Options (defaults = the values used in the benchmark):
%     ell = 0.95, u = 1 (satisfy the compatibility condition for P11-P15),
%     t = 0.95, sigma = 1e-4, rho = 0.5, tauFcn = @(k) 2.^(-k) (k = 0,1,...),
%     initialStep = 'secant', stepMin = 1e-6, stepMax = 1e6,
%     probeEpsilon = 1e-8,
%     plus the protocol options of BM_COMMON_DEFAULTS.

switch variant
    case 'restarted', solverName = 'restarted_hybridplus';
    case 'identity',  solverName = 'hybridplus_cr';
    otherwise
        error('bm:internal', 'Unknown hybridplus variant "%s".', variant);
end

defaults = bm_common_defaults();
defaults.ell = 0.95;
defaults.u = 1;
defaults.t = 0.95;
defaults.sigma = 1e-4;
defaults.rho = 0.5;
defaults.tauFcn = @(k) 2.^(-k);
defaults.initialStep = 'secant';
defaults.stepMin = 1e-6;
defaults.stepMax = 1e6;
defaults.probeEpsilon = 1e-8;
options = bm_merge_options(options, defaults, solverName);
[x, oracle] = bm_validate_common(problem, x0, options, solverName);

bm_require_scalar(options.ell, 'ell', solverName, 0, 1, false, true);
bm_require_scalar(options.u, 'u', solverName, 1, Inf, true, false);
if ~(options.ell < options.u)
    error('bm:invalidOptions', '%s: require ell < u.', solverName);
end
bm_require_scalar(options.t, 't', solverName, 0, 1);
bm_require_scalar(options.sigma, 'sigma', solverName, 0, 1);
bm_require_scalar(options.rho, 'rho', solverName, 0, 1);
if ~isa(options.tauFcn, 'function_handle')
    error('bm:invalidOptions', '%s: tauFcn must be a function handle.', solverName);
end
if ~((ischar(options.initialStep) || (isstring(options.initialStep) && ...
        isscalar(options.initialStep))) && ...
        ismember(char(options.initialStep), {'unit', 'fd', 'secant'}))
    error('bm:invalidOptions', ...
        '%s: initialStep must be ''unit'', ''fd'' or ''secant''.', solverName);
end
initialRule = char(options.initialStep);
bm_require_scalar(options.stepMin, 'stepMin', solverName, 0, 1, false, true);
bm_require_scalar(options.stepMax, 'stepMax', solverName, 1, Inf, true, false);
bm_require_scalar(options.probeEpsilon, 'probeEpsilon', solverName, 0, Inf);

useDiagonal = ~strcmp(variant, 'identity');
t = options.t;
tFactor = sqrt(1 - t^2) / t;

history = bm_history('init', options.maxIterations, options.recordHistory);
k = 0;
restarts = 0;
status = ''; %#ok<NASGU> defensive default; every exit path sets it
message = '';
F = [];

try
    [F, valid] = oracle.evaluate(x);
    if ~valid
        result = bm_make_result(solverName, x, [], 'invalid_initial_point', ...
            'F(x0) is not a finite real vector.', 0, oracle, 0, history, options);
        return;
    end
    history = bm_history('record', history, 0, norm(F), NaN, false);
    xPrev = []; FPrev = []; dPrev = [];

    while true
        if norm(F) <= options.tolerance
            status = 'success';
            break;
        end
        if k >= options.maxIterations
            status = 'max_iterations';
            message = 'iteration limit reached.';
            break;
        end
        tauK = options.tauFcn(k);
        if ~(isscalar(tauK) && isreal(tauK) && isfinite(tauK) && tauK > 0)
            error('bm:invalidOptions', '%s: tauFcn(%d) must be finite and positive.', ...
                solverName, k);
        end

        % ---------------- search direction ----------------
        isRestart = false;
        if k == 0
            direction = -F;                                   % B_0 = I
        else
            s = x - xPrev;
            y = F - FPrev;
            b = ones(size(F));
            if useDiagonal
                active = (s ~= 0);
                b(active) = min(max(y(active) ./ s(active), options.ell), options.u);
            end
            BF = F ./ b;                                      % B_k F_k
            p = dPrev;
            if all(p == 0)
                restartNeeded = true;                         % p = 0
            else
                w = y - ((p' * y) / (p' * p)) * p;
                if all(w == 0)
                    restartNeeded = true;                     % w = 0
                else
                    z = p + (tFactor * norm(p)) * (w / norm(w));
                    numerator = max(0, BF' * z);              % F_k' B_k z
                    denominator = max(p' * z, FPrev' * FPrev);
                    beta = numerator / denominator;
                    direction = -BF + beta * p;
                    restartNeeded = all(direction == 0);      % computed d = 0
                end
            end
            if restartNeeded
                direction = -BF;
                isRestart = true;
                restarts = restarts + 1;
            end
        end
        if any(~isfinite(direction))
            status = 'numerical_error';
            message = 'search direction contains NaN or Inf.';
            break;
        end

        % ---------------- safeguarded initial trial step ----------------
        dd = direction' * direction;
        estimate = 1;
        switch initialRule
            case 'fd'
                epsilon = options.probeEpsilon;
                [Fprobe, validProbe] = oracle.evaluate(x + epsilon * direction);
                if validProbe
                    estimate = abs(F' * direction) / ...
                        abs(direction' * (Fprobe - F) / epsilon);
                end
            case 'secant'
                if k > 0
                    mu = abs(s' * y) / (s' * s);
                    estimate = abs(F' * direction) / (mu * dd);
                end
        end
        if ~isfinite(estimate) || estimate <= 0
            estimate = 1;
        end
        lambda0 = min(max(estimate, options.stepMin), options.stepMax);

        % ---------------- two-sided nonmonotone backtracking ----------------
        fk = 0.5 * (F' * F);
        lambda = lambda0;
        accepted = false;
        while ~accepted
            for sgn = [1, -1]
                [Ft, valid] = oracle.evaluate(x + (sgn * lambda) * direction);
                if valid && 0.5 * (Ft' * Ft) <= (1 + tauK) * fk - options.sigma * lambda^2 * dd
                    accepted = true;
                    xNew = x + (sgn * lambda) * direction;
                    break;
                end
            end
            if ~accepted
                lambda = options.rho * lambda;
                if lambda == 0
                    break;
                end
            end
        end
        if ~accepted
            status = 'line_search_failure';
            message = 'step length underflowed to zero.';
            break;
        end

        xPrev = x; FPrev = F; dPrev = direction;
        x = xNew; F = Ft;
        k = k + 1;
        history = bm_history('record', history, k, norm(F), lambda, isRestart);
    end
catch err
    if strcmp(err.identifier, 'bm:budgetExhausted')
        status = 'max_function_evaluations';
        message = err.message;
    else
        rethrow(err);
    end
end

result = bm_make_result(solverName, x, F, status, message, k, oracle, ...
    restarts, history, options);
end
