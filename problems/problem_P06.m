function F = problem_P06(x)
%PROBLEM_P06 Engval symmetric nonlinear equation (general set).
%   F_1(x) = x_1 (x_1^2 + x_2^2) - 1,
%   F_i(x) = x_i (x_{i-1}^2 + 2 x_i^2 + x_{i+1}^2) - 1,   i = 2, ..., n-1,
%   F_n(x) = x_n (x_{n-1}^2 + x_n^2).
%   The last component has no leading minus sign, as in Li and Fukushima
%   (1999, Problem 2). F is the gradient of (1/4) sum_i (x_i^2+x_{i+1}^2)^2 - sum_{i<n} x_i, so
%   the Jacobian is symmetric positive semidefinite, singular at x = 0.
%   Requires n >= 2.
x = x(:);
n = numel(x);
if n < 2
    error('bm:problemDimension', 'P06 requires n >= 2.');
end
x2 = x.^2;
neighbours = zeros(n, 1);
neighbours(1:n-1) = neighbours(1:n-1) + x2(2:n);
neighbours(2:n) = neighbours(2:n) + x2(1:n-1);
F = x .* (2 * x2 + neighbours) - 1;
F(1) = x(1) * (x2(1) + x2(2)) - 1;
F(n) = x(n) * (x2(n-1) + x2(n));
end
