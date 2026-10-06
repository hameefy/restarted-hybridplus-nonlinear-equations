function F = problem_P14(x)
%PROBLEM_P14 Coupled sine-diffusion problem (SPD set).
%   F(x) = A_4 (x - r) + sin(x) - sin(r),  A_4 = tridiag(-1, 4, -1),
%   r_i = i/(n+1). Exact solution x* = r.
%   J(x) = A_4 + diag(cos x): symmetric, I <= J(x) <= 7 I on R^n.
%   Requires n >= 2.
x = x(:);
n = numel(x);
if n < 2
    error('bm:problemDimension', 'P14 requires n >= 2.');
end
r = (1:n)' / (n + 1);
e = x - r;
F = 4 * e + sin(x) - sin(r);
F(1:n-1) = F(1:n-1) - e(2:n);
F(2:n) = F(2:n) - e(1:n-1);
end
