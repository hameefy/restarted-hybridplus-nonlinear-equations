function F = problem_P04(x)
%PROBLEM_P04 Modified strictly convex function II (general set).
%   F_i(x) = (i/(n+1)) exp(x_i) - 1,   i = 1, ..., n.
x = x(:);
n = numel(x);
F = ((1:n)' / (n + 1)) .* exp(x) - 1;
end
