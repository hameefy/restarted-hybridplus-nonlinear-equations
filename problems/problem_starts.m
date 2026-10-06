function starts = problem_starts(n, seed)
%PROBLEM_STARTS The ten initial points used for every problem.
%
%   STARTS = PROBLEM_STARTS(N, SEED) returns a 1-by-10 struct array with
%   fields id (1..10), name ('guess1', ...) and x0 (N-by-1). The same
%   points are used for every problem and every solver (uniform per class);
%   problem-specific exclusions are applied by the controller using
%   PROBLEM_REGISTRY.excludedStarts.
%
%     guess1  = (1, ..., 1)
%     guess2  = (0.1, ..., 0.1)
%     guess3  = (1 - 1/n, 1 - 2/n, ..., 0)        = ((n-i)/n)_i
%     guess4  = (0, 1/n, 2/n, ..., (n-1)/n)       = ((i-1)/n)_i
%     guess5  = (1, 1/2, ..., 1/n)
%     guess6  = (1 - 2/(n+1), 1 - 4/(n+1), ..., 1 - 2n/(n+1))
%             = ((n + 1 - 2i)/(n + 1))_i          (see below)
%     guess7  = (1/n, 2/n, ..., 1)
%     guess8  = (-1, ..., -1)
%     guess9  =  rand(n,1)   (uniform on (0,1))
%     guess10 = -rand(n,1)   (an independent draw)
%
%   guess6 is a linear ramp from just below +1 to just above -1: it is the only
%   deterministic start with entries of both signs, all its entries lie in
%   (-1, 1) (inside the domain x_i > -1 of P2, so no exclusion is needed),
%   it differs from every other start, and delta'*x0 = -(n-1)/6 ~= 0 for
%   n >= 2, so the P9 formula is defined at it.
%
%   Random starts: a Mersenne-twister stream seeded with SEED + N produces
%   guess9 and then guess10, so the random vectors depend only on (SEED, N)
%   and are identical for all problems and solvers.

if ~(isscalar(n) && isnumeric(n) && n >= 1 && n == floor(n))
    error('bm:invalidDimension', 'n must be a positive integer.');
end
if ~(isscalar(seed) && isnumeric(seed) && seed >= 0 && seed == floor(seed))
    error('bm:invalidSeed', 'seed must be a nonnegative integer.');
end
i = (1:n)';
stream = RandStream('mt19937ar', 'Seed', seed + n);
r9 = rand(stream, n, 1);
r10 = rand(stream, n, 1);
values = {ones(n, 1), 0.1 * ones(n, 1), (n - i) / n, (i - 1) / n, 1 ./ i, ...
    (n + 1 - 2 * i) / (n + 1), i / n, -ones(n, 1), r9, -r10};
starts = struct('id', num2cell(1:10), ...
    'name', arrayfun(@(k) sprintf('guess%d', k), 1:10, 'UniformOutput', false), ...
    'x0', values);
end
