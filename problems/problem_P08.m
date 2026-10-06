function F = problem_P08(x)
%PROBLEM_P08 Modified trigonometric function (general set).
%   F_i(x) = 2 [ n + i (1 - cos x_i) - sin x_i - cos x_i ] (2 sin x_i - cos x_i),
%   i = 1, ..., n. Componentwise (no coupling sum). Level sets are
%   unbounded (periodicity), so the problem lies outside Assumption 1 of
%   the manuscript; it is included as an empirical test.
x = x(:);
n = numel(x);
i = (1:n)';
F = 2 * (n + i .* (1 - cos(x)) - sin(x) - cos(x)) .* (2 * sin(x) - cos(x));
end
