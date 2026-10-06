function varargout = bm_history(action, varargin)
%BM_HISTORY Preallocated per-iteration history used by all solvers.
%
%   H = BM_HISTORY('init', MAXITER, ENABLED) creates a history structure.
%   H = BM_HISTORY('record', H, K, RESNORM, STEP, RESTARTED) stores the
%       residual norm of iterate K (K = 0 is the initial point), the step
%       length that produced it, and whether the direction used to reach
%       it was a restart.
%   H = BM_HISTORY('finalize', H, K) trims the arrays to iterates 0..K.
%
%   Fields: residualNorm (||F(x_k)||), stepLength (lambda_{k-1}, NaN for
%   k = 0), restarted (logical, false for k = 0). When ENABLED is false,
%   only the initial and final residual norms are kept.

switch action
    case 'init'
        maxIter = varargin{1};
        enabled = logical(varargin{2});
        m = 1;
        if enabled, m = maxIter + 1; end
        H = struct('enabled', enabled, ...
            'residualNorm', nan(m, 1), ...
            'stepLength', nan(m, 1), ...
            'restarted', false(m, 1));
        varargout{1} = H;
    case 'record'
        H = varargin{1}; k = varargin{2};
        if H.enabled
            H.residualNorm(k + 1) = varargin{3};
            H.stepLength(k + 1) = varargin{4};
            H.restarted(k + 1) = varargin{5};
        elseif k == 0
            H.residualNorm(1) = varargin{3};
        end
        varargout{1} = H;
    case 'finalize'
        H = varargin{1}; k = varargin{2};
        if H.enabled
            H.residualNorm = H.residualNorm(1:k + 1);
            H.stepLength = H.stepLength(1:k + 1);
            H.restarted = H.restarted(1:k + 1);
        end
        varargout{1} = H;
    otherwise
        error('bm:internal', 'Unknown bm_history action "%s".', action);
end
end
