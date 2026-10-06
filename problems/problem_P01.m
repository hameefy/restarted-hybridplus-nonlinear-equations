function F = problem_P01(x)
%PROBLEM_P01 Modified exponential function (general set).
%   F_1(x) = exp(x_1) - 1,
%   F_i(x) = exp(x_i) + x_i - 1,   i = 2, ..., n.
%   Jacobian: diagonal, symmetric, pointwise positive definite, not
%   uniformly so on R^n.
x = x(:);
F = exp(x) - 1;
F(2:end) = F(2:end) + x(2:end);
end
