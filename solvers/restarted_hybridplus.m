function result = restarted_hybridplus(problem, x0, options)
%RESTARTED_HYBRIDPLUS Restarted hybridplus (Algorithm 1 of the manuscript).
%
%   RESULT = RESTARTED_HYBRIDPLUS(PROBLEM, X0, OPTIONS) solves
%   PROBLEM.F(x) = 0 from X0. See HYBRIDPLUS_CORE for the method, the
%   options and their defaults, and BM_MAKE_RESULT for the result fields.
if nargin < 3, options = struct(); end
result = hybridplus_core(problem, x0, options, 'restarted');
end
