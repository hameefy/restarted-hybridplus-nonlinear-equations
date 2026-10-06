function F = problem_P12(x)
%PROBLEM_P12 Diagonal sine problem (SPD set).
%   F_i(x) = 2 x_i - sin x_i,   i = 1, ..., n.
%   J(x) = diag(2 - cos x): symmetric, I <= J(x) <= 3 I on R^n.
%   Separable with identical components (see manuscript Remark 3.2).
x = x(:);
F = 2 * x - sin(x);
end
