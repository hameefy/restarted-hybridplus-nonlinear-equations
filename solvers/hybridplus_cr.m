function result = hybridplus_cr(problem, x0, options)
%HYBRIDPLUS_CR Restarted hybridplus with B_k = I (diagonal-scaling ablation).
%
%   RESULT = HYBRIDPLUS_CR(PROBLEM, X0, OPTIONS) runs the restarted method
%   with the diagonal matrix replaced by the identity at every iteration,
%   including the restart rule. The options ell and u are accepted for
%   interface uniformity and are not used. See HYBRIDPLUS_CORE.
if nargin < 3, options = struct(); end
result = hybridplus_core(problem, x0, options, 'identity');
end
