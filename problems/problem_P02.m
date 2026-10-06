function F = problem_P02(x)
%PROBLEM_P02 Logarithmic function (general set).
%   F_i(x) = log(x_i + 1) - x_i/n,   i = 1, ..., n,
%   domain D_2 = {x : x_i > -1 for all i}.
%   Outside D_2, MATLAB's log returns -Inf or a complex value; the
%   ResidualOracle then flags the evaluation as invalid and the solver
%   rejects the trial. guess8 = -ones(n,1) is excluded for this problem
%   (it lies on the boundary of D_2).
x = x(:);
n = numel(x);
F = log(x + 1) - x / n;
end
