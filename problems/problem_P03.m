function F = problem_P03(x)
%PROBLEM_P03 Strictly convex function II (general set).
%   F_i(x) = exp(x_i) - 1,   i = 1, ..., n.
%   Separable with identical components: from a constant start every
%   iteration of restarted hybridplus after the first is a restart
%   (manuscript Remark 3.2).
x = x(:);
F = exp(x) - 1;
end
