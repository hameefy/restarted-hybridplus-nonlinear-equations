function result = gnbfgs(problem, x0, options)
%GNBFGS Gauss-Newton-based BFGS method of Li and Fukushima (1999).
%
%   RESULT = GNBFGS(PROBLEM, X0, OPTIONS) implements Algorithm 1 of
%   D.-H. Li and M. Fukushima, SIAM J. Numer. Anal. 37 (1999) 152-172,
%   with g = F:
%     Step 1  stop if ||g_k|| <= tolerance; otherwise solve
%             B_k p + lambda_{k-1}^{-1} (g(x_k + lambda_{k-1} g_k) - g_k) = 0. (12)
%     Step 2  if ||g(x_k + p_k)|| <= rho ||g_k||, take lambda_k = 1.        (13)
%     Step 3  otherwise lambda_k = r^i, i >= 0 the least integer with
%             ||g(x_k+lambda p_k)||^2 - ||g_k||^2 <= -sigma1||lambda g_k||^2
%                          - sigma2||lambda p_k||^2 + eps_k ||g_k||^2.      (11)
%     Step 4  x_{k+1} = x_k + lambda_k p_k.
%     Step 5  s_k = x_{k+1}-x_k, delta_k = g_{k+1}-g_k,
%             y_k = g(x_k + delta_k) - g(x_k); if y_k's_k > 0 apply the BFGS
%             update (14) to B_k, otherwise B_{k+1} = B_k.
%
%   Implementation notes (not changes to the method):
%   * The inverse H_k = B_k^{-1} is stored and updated with the inverse
%     BFGS formula H+ = (I - r s y')H(I - r y s') + r s s', r = 1/(y's),
%     which in exact arithmetic gives the same p_k as solving (12) with the
%     B_k produced by (14). B_0 = H_0 = I (paper, Section 4).
%   * H_k is a dense n-by-n matrix (8 n^2 bytes). Runs with n larger than
%     options.maxDimension are not attempted and return status
%     'not_run_dimension_cap' (counted as failures; see the protocol).
%   * The trial g(x_k + p_k) of Step 2 is reused for i = 0 in Step 3.
%   * If g(x_k + lambda_{k-1} g_k) is not finite/real, p_k cannot be formed
%     and the run stops with status 'invalid_trial_point'. If
%     g(x_k + delta_k) is not finite/real, the update is skipped (as for
%     y's <= 0). Every call of g counts as one residual evaluation.
%
%   Options (defaults = paper values, Section 4):
%     r = 0.1, rho = sqrt(0.9), sigma1 = sigma2 = 1e-5, lambdaMinus1 = 0.01,
%     epsFcn = @(k) 1/(k+1)^2  -- the paper prints eps_k = k^{-2}, which is
%     undefined at k = 0; the index is shifted by one.
%     maxDimension = Inf (BENCHMARK_CONFIG sets the cap).

if nargin < 3, options = struct(); end
solverName = 'gnbfgs';
defaults = bm_common_defaults();
defaults.r = 0.1;
defaults.rho = sqrt(0.9);
defaults.sigma1 = 1e-5;
defaults.sigma2 = 1e-5;
defaults.lambdaMinus1 = 0.01;
defaults.epsFcn = @(k) 1 ./ (k + 1).^2;
defaults.maxDimension = Inf;
options = bm_merge_options(options, defaults, solverName);
[x, oracle] = bm_validate_common(problem, x0, options, solverName);

bm_require_scalar(options.r, 'r', solverName, 0, 1);
bm_require_scalar(options.rho, 'rho', solverName, 0, 1);
bm_require_scalar(options.sigma1, 'sigma1', solverName, 0, Inf);
bm_require_scalar(options.sigma2, 'sigma2', solverName, 0, Inf);
bm_require_scalar(options.lambdaMinus1, 'lambdaMinus1', solverName, 0, Inf);
bm_require_scalar(options.maxDimension, 'maxDimension', solverName, 0, Inf, false, true);
if ~isa(options.epsFcn, 'function_handle')
    error('bm:invalidOptions', '%s: epsFcn must be a function handle.', solverName);
end

n = numel(x);
history = bm_history('init', options.maxIterations, options.recordHistory);
k = 0;
status = ''; %#ok<NASGU> defensive default; every exit path sets it
message = '';
F = [];

if n > options.maxDimension
    result = bm_make_result(solverName, x, [], 'not_run_dimension_cap', ...
        sprintf('n = %d exceeds maxDimension = %g (dense n-by-n matrix).', ...
        n, options.maxDimension), 0, oracle, 0, history, options);
    return;
end

try
    [F, valid] = oracle.evaluate(x);
    if ~valid
        result = bm_make_result(solverName, x, [], 'invalid_initial_point', ...
            'F(x0) is not a finite real vector.', 0, oracle, 0, history, options);
        return;
    end
    history = bm_history('record', history, 0, norm(F), NaN, false);
    H = eye(n);
    lambdaPrev = options.lambdaMinus1;

    while true
        normF = norm(F);
        if normF <= options.tolerance
            status = 'success';
            break;
        end
        if k >= options.maxIterations
            status = 'max_iterations';
            message = 'iteration limit reached.';
            break;
        end
        epsK = options.epsFcn(k);
        if ~(isscalar(epsK) && isreal(epsK) && isfinite(epsK) && epsK > 0)
            error('bm:invalidOptions', '%s: epsFcn(%d) must be finite and positive.', ...
                solverName, k);
        end

        % ---- Step 1: p_k = -H_k * lambda_{k-1}^{-1}(g(x_k + lambda_{k-1} g_k) - g_k)
        [Fshift, valid] = oracle.evaluate(x + lambdaPrev * F);
        if ~valid
            status = 'invalid_trial_point';
            message = 'g(x_k + lambda_{k-1} g_k) is not finite/real; p_k undefined.';
            break;
        end
        v = (Fshift - F) / lambdaPrev;
        p = -(H * v);
        if any(~isfinite(p))
            status = 'numerical_error';
            message = 'search direction contains NaN or Inf.';
            break;
        end

        % ---- Step 2: unit step test (13)
        [Funit, validUnit] = oracle.evaluate(x + p);
        if validUnit && norm(Funit) <= options.rho * normF
            lambda = 1;
            Fnew = Funit;
        else
            % ---- Step 3: backtracking on (11), i = 0 reuses g(x_k + p_k)
            gg = normF^2;
            pp = p' * p;
            lambda = 1;
            accepted = false;
            while ~accepted
                if lambda == 1
                    Ft = Funit; valid = validUnit;
                else
                    [Ft, valid] = oracle.evaluate(x + lambda * p);
                end
                if valid && (Ft' * Ft) - gg <= -options.sigma1 * lambda^2 * gg ...
                        - options.sigma2 * lambda^2 * pp + epsK * gg
                    accepted = true;
                    Fnew = Ft;
                else
                    lambda = options.r * lambda;
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
        end

        % ---- Steps 4-5: new iterate and BFGS update
        s = lambda * p;
        xNew = x + s;
        delta = Fnew - F;
        [Fdelta, validDelta] = oracle.evaluate(x + delta);
        if validDelta
            y = Fdelta - F;
            ys = y' * s;
            if ys > 0
                Hy = H * y;
                rr = 1 / ys;
                u1 = (0.5 * (rr^2 * (y' * Hy) + rr)) * s - rr * Hy;
                H = H + (u1 * s' + s * u1');
            end
        end

        x = xNew; F = Fnew; lambdaPrev = lambda;
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
