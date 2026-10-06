function defaults = bm_common_defaults()
%BM_COMMON_DEFAULTS Protocol options shared by every solver.
%
%   tolerance              success iff ||F(x_k)||_2 <= tolerance
%   maxIterations          maximum number of accepted iterations
%   maxFunctionEvaluations maximum number of residual evaluations (every
%                          call of F counts, including the initial one,
%                          rejected trials and derivative-free probes)
%   recordHistory          store ||F(x_k)|| and step data per iteration
%
%   These values define the stopping rule of the benchmark: a run is a
%   failure when ||F(x_k)|| > 1e-6 after 1000 iterations or
%   2000 residual evaluations.
defaults = struct( ...
    'tolerance', 1e-6, ...
    'maxIterations', 1000, ...
    'maxFunctionEvaluations', 2000, ...
    'recordHistory', true);
end
