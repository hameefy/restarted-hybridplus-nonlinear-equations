function F = problem_P09(x)
%PROBLEM_P09 Function derived from the Chandrasekhar H-equation (general set).
%   For the H-equation see Kelley (1980). The function computed is
%       F_i(x) = x_i - 4 / ( c n (delta'x) sum_j (1/delta_j) ),  i = 1..n,
%   with delta_j = (j - 0.5)/n and c = 0.9, as stated in the manuscript.
%   In the code below, with x a column vector, the expression u(i)*x(j) is
%   the inner product delta'*x (a scalar) and c/2*n evaluates to c*n/2, so
%       A    = 1 - (c*n/2) * (delta'*x) * sum_j 1/(2 delta_j),
%       F(x) = x - (1 - A)^(-1) * ones(n,1),
%   which is the formula above. It is not the componentwise H-equation.
%   Cost O(n) time and memory. F is undefined where delta'*x = 0; there
%   the value is Inf and the oracle rejects the trial.
x = x(:);                 % the code requires a column vector
n = length(x);
F = zeros(n, 1);
i = 1:(n);
c = 0.9;
u = zeros(1, n);          % preallocation only
u(i) = (i - 0.5) / n;
j = 1:(n);
A = 1. - ((c/2*n) * sum((u(i) * x(j)) ./ (u(i) + u(j))));
F(i) = x(i) - (1. - A).^(-1);
F = F(:);
end
