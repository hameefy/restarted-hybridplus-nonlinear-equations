function F = problem_P07(x)
%PROBLEM_P07 Trigexp function (general set).
%   F_1(x) = 3 x_1^3 + 2 x_2 - 5 + sin(x_1 - x_2) sin(x_1 + x_2),
%   F_i(x) = -x_{i-1} exp(x_{i-1} - x_i) + x_i (4 + 3 x_i^2) + 2 x_{i+1}
%            + sin(x_i - x_{i+1}) sin(x_i + x_{i+1}) - 8,   i = 2, ..., n-1,
%   F_n(x) = -x_{n-1} exp(x_{n-1} - x_n) + 4 x_n - 3.
%   Requires n >= 2. Jacobian generally nonsymmetric.
x = x(:);
n = numel(x);
if n < 2
    error('bm:problemDimension', 'P07 requires n >= 2.');
end
F = zeros(n, 1);
F(1) = 3 * x(1)^3 + 2 * x(2) - 5 + sin(x(1) - x(2)) * sin(x(1) + x(2));
if n > 2
    im = 1:n-2;           % i-1
    ii = 2:n-1;           % i
    ip = 3:n;             % i+1
    F(ii) = -x(im) .* exp(x(im) - x(ii)) + x(ii) .* (4 + 3 * x(ii).^2) ...
        + 2 * x(ip) + sin(x(ii) - x(ip)) .* sin(x(ii) + x(ip)) - 8;
end
F(n) = -x(n-1) * exp(x(n-1) - x(n)) + 4 * x(n) - 3;
end
