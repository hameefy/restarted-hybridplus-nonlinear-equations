function [x0, oracle] = bm_validate_common(problem, x0, options, solverName)
%BM_VALIDATE_COMMON Check the common solver inputs and build the oracle.
%
%   [X0, ORACLE] = BM_VALIDATE_COMMON(PROBLEM, X0, OPTIONS, SOLVERNAME)
%   verifies that PROBLEM is a structure with a function handle field F,
%   that X0 is a finite real vector, and that the protocol options in
%   OPTIONS are admissible. It returns X0 as a column vector and a
%   ResidualOracle with the evaluation budget of OPTIONS.

if ~isstruct(problem) || ~isscalar(problem) || ~isfield(problem, 'F') || ...
        ~isa(problem.F, 'function_handle')
    error('bm:invalidProblem', '%s: problem.F must be a function handle.', solverName);
end
if ~isnumeric(x0) || ~isreal(x0) || isempty(x0) || ~isvector(x0) || ...
        any(~isfinite(x0(:)))
    error('bm:invalidStart', '%s: x0 must be a nonempty finite real vector.', solverName);
end
x0 = double(x0(:));

bm_require_scalar(options.tolerance, 'tolerance', solverName, 0, Inf, true);
bm_require_integer(options.maxIterations, 'maxIterations', solverName);
bm_require_integer(options.maxFunctionEvaluations, 'maxFunctionEvaluations', solverName);
if ~(islogical(options.recordHistory) || isnumeric(options.recordHistory)) || ...
        ~isscalar(options.recordHistory)
    error('bm:invalidOptions', '%s: recordHistory must be a logical scalar.', solverName);
end

oracle = ResidualOracle(problem.F, numel(x0), options.maxFunctionEvaluations);
end

function bm_require_integer(value, name, solverName)
if ~(isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value) && ...
        value >= 1 && value == floor(value))
    error('bm:invalidOptions', '%s: %s must be a positive integer.', solverName, name);
end
end
