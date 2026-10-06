function F = problem_P10(x)
%PROBLEM_P10 Broyden tridiagonal function (general set).
%   F_1(x) = (3 - h x_1) x_1 - 2 x_2 + 1,
%   F_i(x) = (3 - h x_i) x_i - x_{i-1} - 2 x_{i+1} + 1,   i = 2, ..., n-1,
%   F_n(x) = (3 - h x_n) x_n - x_{n-1} + 1,
%   with the standard coefficient h = 2.
%   Requires n >= 2. Jacobian nonsymmetric (off-diagonals -2 and -1).
x = x(:);
n = numel(x);
if n < 2
    error('bm:problemDimension', 'P10 requires n >= 2.');
end
h = 2;
F = (3 - h * x) .* x + 1;
F(1:n-1) = F(1:n-1) - 2 * x(2:n);
F(2:n) = F(2:n) - x(1:n-1);
end
