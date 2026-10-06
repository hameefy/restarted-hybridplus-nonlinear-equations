function F = problem_P13(x)
%PROBLEM_P13 Symmetric linear tridiagonal control problem (SPD set).
%   F_1(x) = (5/2) x_1 + x_2 - 1,
%   F_i(x) = x_{i-1} + (5/2) x_i + x_{i+1} - 1,   i = 2, ..., n-1,
%   F_n(x) = x_{n-1} + (5/2) x_n - 1.
%   J = tridiag(1, 5/2, 1): symmetric, (1/2) I <= J <= (9/2) I.
%   Requires n >= 2.
x = x(:);
n = numel(x);
if n < 2
    error('bm:problemDimension', 'P13 requires n >= 2.');
end
F = 2.5 * x - 1;
F(1:n-1) = F(1:n-1) + x(2:n);
F(2:n) = F(2:n) + x(1:n-1);
end
