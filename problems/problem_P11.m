function F = problem_P11(x)
%PROBLEM_P11 Symmetric sine boundary-value problem (SPD set).
%   F(x) = A_8 x + sin(x) - 1,  A_8 = tridiag(-1, 8, -1)  (n-by-n).
%   J(x) = A_8 + diag(cos x): symmetric, 5 I <= J(x) <= 11 I on R^n.
%   Requires n >= 2.
x = x(:);
n = numel(x);
if n < 2
    error('bm:problemDimension', 'P11 requires n >= 2.');
end
F = 8 * x + sin(x) - 1;
F(1:n-1) = F(1:n-1) - x(2:n);
F(2:n) = F(2:n) - x(1:n-1);
end
