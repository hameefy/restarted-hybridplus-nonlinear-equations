function F = problem_P15(x)
%PROBLEM_P15 Coupled hyperbolic-tangent diffusion problem (SPD set).
%   F(x) = A_4 (x - r) + tanh(x) - tanh(r),  A_4 = tridiag(-1, 4, -1),
%   r_i = i/(n+1). Exact solution x* = r.
%   J(x) = A_4 + diag(sech^2 x): symmetric, 2 I <= J(x) <= 7 I on R^n.
%   Requires n >= 2.
x = x(:);
n = numel(x);
if n < 2
    error('bm:problemDimension', 'P15 requires n >= 2.');
end
r = (1:n)' / (n + 1);
e = x - r;
F = 4 * e + tanh(x) - tanh(r);
F(1:n-1) = F(1:n-1) - e(2:n);
F(2:n) = F(2:n) - e(1:n-1);
end
