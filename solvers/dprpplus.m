function result = dprpplus(problem, x0, options)
%DPRPPLUS Diagonal PRP+ method of Mohammad and Awwal (2022).
%
%   RESULT = DPRPPLUS(PROBLEM, X0, OPTIONS) implements Algorithm 1 (dprp+)
%   of H. Mohammad and A. M. Awwal, Numer. Algorithms 91 (2022) 1441-1460:
%     d_0 = -F_0,  d_k = -Q_k F_k + beta_k d_{k-1},
%     q_k^i = max{ min{ y_i/s_i, u }, ell } if s_i ~= 0, 1 otherwise   (9)
%     Q_k = diag(1./q_k),
%     beta_k = max{0, <F_k, y_{k-1}>} / max{ ||F_{k-1}||^2, ||d_{k-1}|| }  (10)
%   and the one-sided nonmonotone line search (12): lambda = theta*rho^i,
%   i the least nonnegative integer with
%     f(x_k + lambda d_k) <= (1 + tau_k) f(x_k) - sigma lambda^2 ||d_k||^2.
%
%   NOTE on (10): the denominator is transcribed exactly as printed in the
%   paper, i.e. with the NORM ||d_{k-1}|| (not d_{k-1}'y_{k-1}).
%
%   Options (defaults = paper values for theta, rho; values chosen by the
%   authors where the paper does not state them):
%     theta = 1, rho = 0.2                       (paper, Section 4)
%     sigma = 0.01                               (not stated in the paper)
%     tauFcn = @(k) exp(-(k+1).^(k+1))           (our reading of the
%                                                 printed 1/e^{k+1^{k+1}})
%     ell = 1e-10, u = 1e10                      (not stated in the paper)
%   In double precision tauFcn(k) = 0 for k >= 4 (e^{-3125} underflows);
%   the line search is then monotone. This is a property of the prescribed
%   sequence, not an implementation choice.

if nargin < 3, options = struct(); end
solverName = 'dprpplus';
defaults = bm_common_defaults();
defaults.ell = 1e-10;
defaults.u = 1e10;
defaults.theta = 1;
defaults.rho = 0.2;
defaults.sigma = 0.01;
defaults.tauFcn = @(k) exp(-(k + 1).^(k + 1));
options = bm_merge_options(options, defaults, solverName);
[x, oracle] = bm_validate_common(problem, x0, options, solverName);

bm_require_scalar(options.ell, 'ell', solverName, 0, 1, false, true);
bm_require_scalar(options.u, 'u', solverName, 1, Inf, true, false);
bm_require_scalar(options.theta, 'theta', solverName, 0, 1, false, true);
bm_require_scalar(options.rho, 'rho', solverName, 0, 1);
bm_require_scalar(options.sigma, 'sigma', solverName, 0, 1);
if ~isa(options.tauFcn, 'function_handle')
    error('bm:invalidOptions', '%s: tauFcn must be a function handle.', solverName);
end

history = bm_history('init', options.maxIterations, options.recordHistory);
k = 0;
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
        if ~(isscalar(tauK) && isreal(tauK) && isfinite(tauK) && tauK >= 0)
            error('bm:invalidOptions', '%s: tauFcn(%d) must be finite and nonnegative.', ...
                solverName, k);
        end

        if k == 0
            direction = -F;                                   % Q_0 = I
        else
            s = x - xPrev;
            y = F - FPrev;
            q = ones(size(F));
            active = (s ~= 0);
            q(active) = max(min(y(active) ./ s(active), options.u), options.ell);
            denominator = max(FPrev' * FPrev, norm(dPrev));   % eq. (10) as printed
            beta = max(0, F' * y) / denominator;
            direction = -(F ./ q) + beta * dPrev;
        end
        if any(~isfinite(direction))
            status = 'numerical_error';
            message = 'search direction contains NaN or Inf.';
            break;
        end

        fk = 0.5 * (F' * F);
        dd = direction' * direction;
        lambda = options.theta;
        accepted = false;
        while ~accepted
            [Ft, valid] = oracle.evaluate(x + lambda * direction);
            if valid && 0.5 * (Ft' * Ft) <= (1 + tauK) * fk - options.sigma * lambda^2 * dd
                accepted = true;
            else
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
        x = x + lambda * direction; F = Ft;
        k = k + 1;
        history = bm_history('record', history, k, norm(F), lambda, false);
    end
catch err
    if strcmp(err.identifier, 'bm:budgetExhausted')
        status = 'max_function_evaluations';
        message = err.message;
    else
        rethrow(err);
    end
end

result = bm_make_result(solverName, x, F, status, message, k, oracle, 0, ...
    history, options);
end
