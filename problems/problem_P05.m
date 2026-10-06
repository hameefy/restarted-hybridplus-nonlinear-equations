function F = problem_P05(x)
%PROBLEM_P05 Tridiagonal exponential function (general set).
%   With h = 1/(n+1):
%   F_1(x) = x_1 - exp(cos(h (x_1 + x_2))),
%   F_i(x) = x_i - exp(cos(h (x_{i-1} + x_i + x_{i+1}))),  i = 2, ..., n-1,
%   F_n(x) = x_n - exp(cos(h (x_{n-1} + x_n))).
%   Requires n >= 2. Jacobian generally nonsymmetric.
x = x(:);
n = numel(x);
if n < 2
    error('bm:problemDimension', 'P05 requires n >= 2.');
end
h = 1 / (n + 1);
sums = x;
sums(1:n-1) = sums(1:n-1) + x(2:n);
sums(2:n) = sums(2:n) + x(1:n-1);
F = x - exp(cos(h * sums));
end
