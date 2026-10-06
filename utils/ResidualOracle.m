classdef ResidualOracle < handle
%RESIDUALORACLE Counted, budget-enforcing wrapper around a residual map F.
%
%   ORACLE = RESIDUALORACLE(FUN, N, BUDGET) wraps the function handle FUN,
%   which must map a column vector x in R^N to a vector F(x) in R^N.
%
%   [F, VALID] = ORACLE.evaluate(X) evaluates FUN(X) and increments the
%   evaluation counter. Every call counts, including rejected line-search
%   trials and derivative-free probes. VALID is false when FUN(X) is not a
%   finite real vector of length N (for example, outside the domain of a
%   logarithm); the solver must then treat the trial as rejected. F is
%   returned as a column vector (empty when VALID is false).
%
%   When the counter has already reached BUDGET, evaluate() does NOT call
%   FUN and throws an error with identifier 'bm:budgetExhausted'. Solvers
%   catch this identifier and report status 'max_function_evaluations'.
%   Exceptions thrown by FUN itself are not caught here: they indicate an
%   implementation error and propagate to the benchmark controller.
%
%   This class is shared by every solver so that residual evaluations are
%   counted identically for all methods.

    properties (SetAccess = private)
        Fun            % function handle x -> F(x)
        Dimension      % problem dimension n
        Budget         % maximum number of evaluations
        Count = 0      % evaluations performed so far
        InvalidCount = 0  % evaluations that returned non-finite/complex values
    end

    methods
        function obj = ResidualOracle(fun, n, budget)
            if ~isa(fun, 'function_handle')
                error('bm:invalidOracle', 'FUN must be a function handle.');
            end
            if ~(isscalar(n) && isnumeric(n) && n >= 1 && n == floor(n))
                error('bm:invalidOracle', 'N must be a positive integer.');
            end
            if ~(isscalar(budget) && isnumeric(budget) && budget >= 1 && ...
                    budget == floor(budget))
                error('bm:invalidOracle', 'BUDGET must be a positive integer.');
            end
            obj.Fun = fun;
            obj.Dimension = n;
            obj.Budget = budget;
        end

        function [F, valid] = evaluate(obj, x)
            if obj.Count >= obj.Budget
                error('bm:budgetExhausted', ...
                    'Residual-evaluation budget of %d reached.', obj.Budget);
            end
            obj.Count = obj.Count + 1;
            F = obj.Fun(x);
            valid = isnumeric(F) && isreal(F) && numel(F) == obj.Dimension ...
                && all(isfinite(F(:)));
            if valid
                F = double(F(:));
            else
                F = [];
                obj.InvalidCount = obj.InvalidCount + 1;
            end
        end

        function r = remaining(obj)
            r = obj.Budget - obj.Count;
        end
    end
end
