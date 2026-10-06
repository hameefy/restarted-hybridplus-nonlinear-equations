function result = bm_make_result(solverName, x, F, status, message, ...
    iterations, oracle, restarts, history, options)
%BM_MAKE_RESULT Assemble the uniform solver result structure.
%
%   Fields of RESULT:
%     solver               solver name
%     x                    final iterate (column vector)
%     finalResidual        ||F(x)||_2 (Inf if F is unavailable)
%     status               termination status (see README, failure taxonomy)
%     success              true iff finalResidual <= options.tolerance
%     message              human-readable termination reason
%     iterations           number of accepted iterations
%     functionEvaluations  residual evaluations (every call of F)
%     invalidEvaluations   evaluations that returned non-finite/complex values
%     jacobianEvaluations  always 0 (all methods are derivative free)
%     restarts             number of restarts (hybridplus variants; else 0)
%     history              per-iteration history (see BM_HISTORY)
%
%   Success is decided ONLY by the residual test, consistently with the
%   protocol: a run whose final residual satisfies the tolerance is a
%   success, every other run is a failure whatever its status.

if isempty(F)
    finalResidual = Inf;
else
    finalResidual = norm(F);
end
success = isfinite(finalResidual) && finalResidual <= options.tolerance;
if success
    status = 'success';
    message = '';
end
result = struct( ...
    'solver', solverName, ...
    'x', x, ...
    'finalResidual', finalResidual, ...
    'status', status, ...
    'success', success, ...
    'message', message, ...
    'iterations', iterations, ...
    'functionEvaluations', oracle.Count, ...
    'invalidEvaluations', oracle.InvalidCount, ...
    'jacobianEvaluations', 0, ...
    'restarts', restarts, ...
    'history', bm_history('finalize', history, iterations));
end
