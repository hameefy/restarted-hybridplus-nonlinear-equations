function passed = run_self_checks()
%RUN_SELF_CHECKS Correctness checks of the framework (seconds, n <= 20).
%
%   PASSED = RUN_SELF_CHECKS() verifies the problems against independent
%   formulas, the restart logic on an exact scalar example (F = 1.5x, in
%   which w = 0 at k = 1), the evaluation budget, option validation,
%   the performance-profile computation on a hand-computed case, and runs
%   the complete pipeline (controller, tables, profiles, plots) on a tiny
%   configuration in a temporary folder that is deleted afterwards.
%   These are software tests; they produce no benchmark results.

setup_paths();
checks = {
    @check_problems
    @check_scalar_example
    @check_oracle_budget
    @check_option_validation
    @check_solvers_small
    @check_dimension_cap
    @check_initial_step_rules
    @check_profile_computation
    @check_pipeline
    };
passed = true;
for i = 1:numel(checks)
    name = func2str(checks{i});
    try
        checks{i}();
        fprintf('[PASS] %s\n', name);
    catch err
        passed = false;
        fprintf(2, '[FAIL] %s: %s\n', name, err.message);
    end
end
if passed
    fprintf('All self checks passed.\n');
else
    fprintf(2, 'Some self checks FAILED; do not start the benchmark.\n');
end
end

% ======================================================================
function check_problems()
registry = problem_registry();
n = 7;
x = linspace(0.2, 0.9, n)';
for e = registry
    F = e.F(x);
    assert(isequal(size(F), [n 1]) && all(isfinite(F)) && isreal(F), ...
        '%s does not return a finite real %d-by-1 vector.', e.name, n);
end
e1 = ones(n - 1, 1);
tri = @(a, b, c) diag(b * ones(n, 1)) + diag(a * e1, -1) + diag(c * e1, 1);
% P01
ref = exp(x) + x - 1; ref(1) = exp(x(1)) - 1;
assert(norm(problem_P01(x) - ref) < 1e-14, 'P01 mismatch.');
% P05
h = 1 / (n + 1); S = tri(1, 1, 1) * x;
assert(norm(problem_P05(x) - (x - exp(cos(h * S)))) < 1e-14, 'P05 mismatch.');
% P06 = gradient of Psi(x) = 1/4 sum (x_i^2+x_{i+1}^2)^2 - sum_{i<n} x_i
psi = @(v) 0.25 * sum((v(1:end-1).^2 + v(2:end).^2).^2) - sum(v(1:end-1));
g = zeros(n, 1); step = 1e-6;
for i = 1:n
    ei = zeros(n, 1); ei(i) = step;
    g(i) = (psi(x + ei) - psi(x - ei)) / (2 * step);
end
assert(norm(problem_P06(x) - g) < 1e-7, 'P06 is not the gradient of the Engval function.');
% P07 (componentwise reference)
ref = zeros(n, 1);
ref(1) = 3*x(1)^3 + 2*x(2) - 5 + sin(x(1)-x(2))*sin(x(1)+x(2));
for i = 2:n-1
    ref(i) = -x(i-1)*exp(x(i-1)-x(i)) + x(i)*(4+3*x(i)^2) + 2*x(i+1) ...
        + sin(x(i)-x(i+1))*sin(x(i)+x(i+1)) - 8;
end
ref(n) = -x(n-1)*exp(x(n-1)-x(n)) + 4*x(n) - 3;
assert(norm(problem_P07(x) - ref) < 1e-12, 'P07 mismatch.');
% P09 closed form (see problem_P09.m)
delta = ((1:n)' - 0.5) / n; c = 0.9;
ref = x - 4 / (c * n * (delta' * x) * sum(1 ./ delta));
assert(norm(problem_P09(x) - ref) < 1e-12, 'P09 differs from its documented closed form.');
% P10
ref = (3 - 2 * x) .* x + 1 - tri(1, 0, 2) * x;
assert(norm(problem_P10(x) - ref) < 1e-14, 'P10 mismatch.');
% P11, P13 (matrix references)
assert(norm(problem_P11(x) - (tri(-1, 8, -1) * x + sin(x) - 1)) < 1e-14, 'P11 mismatch.');
assert(norm(problem_P13(x) - (tri(1, 2.5, 1) * x - 1)) < 1e-14, 'P13 mismatch.');
% P14, P15 vanish at r
r = (1:n)' / (n + 1);
assert(norm(problem_P14(r)) < 1e-14 && norm(problem_P15(r)) < 1e-14, ...
    'P14/P15 do not vanish at their known solution.');
% SPD bounds of P11-P15 at a random point
rs = RandStream('mt19937ar', 'Seed', 1);
z = 4 * randn(rs, n, 1);
J11 = tri(-1, 8, -1) + diag(cos(z));
J13 = tri(1, 2.5, 1);
J14 = tri(-1, 4, -1) + diag(cos(z));
J15 = tri(-1, 4, -1) + diag(sech(z).^2);
bounds = {J11, 5, 11; diag(2 - cos(z)), 1, 3; J13, 0.5, 4.5; J14, 1, 7; J15, 2, 7};
for i = 1:size(bounds, 1)
    ev = eig(bounds{i, 1});
    assert(min(ev) >= bounds{i, 2} - 1e-12 && max(ev) <= bounds{i, 3} + 1e-12, ...
        'Spectral bounds of P%d violated.', 10 + i);
end
% start points
st = problem_starts(n, 20260930);
assert(numel(st) == 10 && all(st(9).x0 > 0) && all(st(10).x0 < 0), ...
    'Start points malformed.');
assert(norm(st(6).x0 - (1 - 2 * (1:n)' / (n + 1))) < 1e-15 && ...
    all(abs(st(6).x0) < 1), 'guess6 differs from (n+1-2i)/(n+1) or leaves (-1,1).');
X = [st.x0];
for a = 1:9
    for b = a+1:10
        assert(~isequal(X(:, a), X(:, b)), 'guess%d and guess%d coincide.', a, b);
    end
end
assert(abs(delta' * st(6).x0 + (n - 1) / 6) < 1e-12, ...
    'delta''*guess6 differs from -(n-1)/6 (P9 would be undefined if it were 0).');
end

function check_scalar_example()
% F(x) = 1.5x, x0 = 1, unit initial trial: x1 = -0.5 (negative trial
% accepted); at k = 1, y is parallel to p, so w = 0 and the restart gives
% d1 = -B1 F1 = 0.5 with B1 = 1/1.5, hence x2 = 0 exactly.
problem = struct('F', @(x) 1.5 * x);
opts = struct('tolerance', 0, 'ell', 1, 'u', 1.5, 'rho', 0.5, 'sigma', 0.1, ...
    'tauFcn', @(k) 0.1 * 2.^(-k), 't', 0.5, 'initialStep', 'unit');
res = restarted_hybridplus(problem, 1, opts);
assert(res.success && res.iterations == 2 && res.restarts == 1 && res.x == 0, ...
    'Restarted hybridplus does not reproduce the scalar example (x1=-0.5, restart, x2=0).');
assert(isequal(res.history.residualNorm, [1.5; 0.75; 0]) && ...
    isequal(res.history.restarted, [false; false; true]), 'Scalar example history mismatch.');
end

function check_oracle_budget()
oracle = ResidualOracle(@(x) x, 2, 2);
oracle.evaluate([1; 2]);
[~, valid] = oracle.evaluate([NaN; 2]);
assert(~valid && oracle.InvalidCount == 1, 'Invalid values not flagged.');
threw = false;
try
    oracle.evaluate([1; 2]);
catch err
    threw = strcmp(err.identifier, 'bm:budgetExhausted');
end
assert(threw && oracle.Count == 2, 'Budget not enforced.');
problem = struct('F', @problem_P05);
names = {'restarted_hybridplus', 'hybridplus_cr', 'dprpplus', 'gnbfgs'};
for i = 1:numel(names)
    res = feval(names{i}, problem, 3 * ones(20, 1), struct('maxFunctionEvaluations', 4));
    assert(res.functionEvaluations <= 4, '%s exceeded the evaluation budget.', names{i});
    assert(res.success || strcmp(res.status, 'max_function_evaluations'), ...
        '%s: unexpected status %s under a tiny budget.', names{i}, res.status);
end
end

function check_option_validation()
problem = struct('F', @problem_P13);
threw = false;
try
    restarted_hybridplus(problem, ones(5, 1), struct('sigmaa', 0.1));
catch err
    threw = strcmp(err.identifier, 'bm:invalidOptions');
end
assert(threw, 'A misspelled option was not rejected.');
threw = false;
try
    restarted_hybridplus(problem, ones(5, 1), struct('t', 1));
catch err
    threw = strcmp(err.identifier, 'bm:invalidOptions');
end
assert(threw, 't = 1 was not rejected.');
cfg = benchmark_config();
assert(strcmp(cfg.initialStepRule, 'secant') && cfg.t == 0.95, ...
    'benchmark_config does not use the secant rule and t = 0.95.');
names = {'restarted_hybridplus', 'hybridplus_cr'};
for i = 1:numel(names)
    o = cfg.solvers(strcmp({cfg.solvers.name}, names{i})).options;
    assert(o.t == 0.95 && strcmp(o.initialStep, 'secant'), ...
        '%s is not configured with t = 0.95 and the secant rule.', names{i});
end
assert(isequal(sort({cfg.solvers([cfg.solvers.enabled]).name}), ...
    sort({'restarted_hybridplus', 'hybridplus_cr', 'dprpplus', 'gnbfgs'})), ...
    'The enabled solvers are not the four solvers of the benchmark.');
bm_validate_config(cfg);
bad = cfg; bad.solvers(1).options.tolerance = 1e-3;
threw = false;
try
    bm_validate_config(bad);
catch err
    threw = strcmp(err.identifier, 'bm:config');
end
assert(threw, 'A solver-specific protocol override was not rejected.');
end

function check_solvers_small()
cfg = benchmark_config();
problem = struct('F', @problem_P13);
for s = 1:numel(cfg.solvers)
    opts = cfg.solvers(s).options;
    res = cfg.solvers(s).fun(problem, ones(10, 1), opts);
    validStatus = {'success', 'max_iterations', 'max_function_evaluations', ...
        'line_search_failure', 'numerical_error', 'invalid_trial_point'};
    assert(ismember(res.status, validStatus), '%s returned status %s.', ...
        cfg.solvers(s).name, res.status);
    assert(res.success == (res.finalResidual <= 1e-6), ...
        '%s: success flag inconsistent with the residual test.', cfg.solvers(s).name);
    assert(res.functionEvaluations <= 2000 && res.iterations <= 1000, ...
        '%s violates the protocol limits.', cfg.solvers(s).name);
    fprintf('       %-22s P13 n=10: %-26s it=%d fe=%d res=%.2e\n', ...
        cfg.solvers(s).name, res.status, res.iterations, ...
        res.functionEvaluations, res.finalResidual);
    assert(numel(res.history.residualNorm) == res.iterations + 1, ...
        '%s history length mismatch.', cfg.solvers(s).name);
end
end

function check_dimension_cap()
res = gnbfgs(struct('F', @problem_P13), ones(20, 1), struct('maxDimension', 10));
assert(strcmp(res.status, 'not_run_dimension_cap') && res.functionEvaluations == 0 ...
    && ~res.success, 'gnbfgs dimension cap not honoured.');
end

function check_initial_step_rules()
problem = struct('F', @problem_P13);
x0 = ones(10, 1);
% the defaults must be the benchmark values (secant rule, t = 0.95)
a = restarted_hybridplus(problem, x0, struct());
b = restarted_hybridplus(problem, x0, struct('initialStep', 'secant', 't', 0.95));
assert(isequal(a.iterations, b.iterations) && isequal(a.functionEvaluations, ...
    b.functionEvaluations) && isequal(a.x, b.x), 'Defaults differ from secant, t = 0.95.');
a = restarted_hybridplus(problem, x0, struct('initialStep', 'unit'));
rules = {'unit', 'secant', 'fd'};
for i = 1:numel(rules)
    res = restarted_hybridplus(problem, x0, struct('initialStep', rules{i}));
    assert(res.functionEvaluations <= 2000 && res.iterations <= 1000 && ...
        res.success == (res.finalResidual <= 1e-6), 'Rule %s violates the protocol.', rules{i});
    steps = res.history.stepLength(2:end);
    assert(all(steps > 0 & steps <= 1e6), 'Rule %s: accepted step outside (0, stepMax].', rules{i});
    fprintf('       initialStep=%-7s P13 n=10: %-26s it=%d fe=%d\n', rules{i}, ...
        res.status, res.iterations, res.functionEvaluations);
end
% 'fd' spends one probe per iteration: at least 2 evaluations per iteration
res = restarted_hybridplus(problem, x0, struct('initialStep', 'fd'));
assert(res.functionEvaluations >= 2 * res.iterations + 1, 'Rule fd: probes are not counted.');
% clipping: with stepMin = stepMax = 1 every rule equals 'unit'
c = restarted_hybridplus(problem, x0, struct('initialStep', 'secant', ...
    'stepMin', 1, 'stepMax', 1));
assert(isequal(c.functionEvaluations, a.functionEvaluations) && isequal(c.x, a.x), ...
    'Clipping to [1,1] does not reproduce the unit rule.');
threw = false;
try
    restarted_hybridplus(problem, x0, struct('initialStep', 'bb'));
catch err
    threw = strcmp(err.identifier, 'bm:invalidOptions');
end
assert(threw, 'An unknown initialStep rule was not rejected.');
end

function check_profile_computation()
% Hand-computed example: 3 cases, 2 solvers.
%   case A: s1 = 10, s2 = 20 -> ratios 1, 2
%   case B: s1 fails, s2 = 5 -> Inf, 1
%   case C: both fail         -> Inf, Inf
T = table(["A"; "A"; "B"; "B"; "C"; "C"], ["s1"; "s2"; "s1"; "s2"; "s1"; "s2"], ...
    logical([1; 1; 0; 1; 0; 0]), [10; 20; 3; 5; 7; 7], ...
    'VariableNames', {'case_id', 'solver', 'success', 'iterations'});
p = compute_performance_profile(T, 'iterations', {'s1', 's2'}, 1);
assert(isequal(p.tau, [1; 2]), 'Profile breakpoints wrong.');
assert(isequal(p.rho, [1/3, 1/3; 1/3, 2/3]), 'Profile values wrong.');
assert(isequal(p.solvedFraction, [1/3, 2/3]), 'Solved fractions wrong.');
end

function check_pipeline()
cfg = benchmark_config();
cfg.startIds = [1 8 9];
cfg.solvers(strcmp({cfg.solvers.name}, 'gnbfgs')).options.maxDimension = 11;
cfg.dimensions = [10 12];            % gnbfgs runs at n = 10 only
cfg.verbose = false;
cfg.analysis.historyCases = {'P13_n10_guess1'};
cfg.outputDirectory = fullfile(tempdir, sprintf('rhp_selfcheck_%s', ...
    char(datetime('now', 'Format', 'yyyyMMddHHmmssSSS'))));
cleanup = onCleanup(@() rmdir(cfg.outputDirectory, 's'));
write_protocol(cfg);
rawFile = benchmark_controller(cfg);
T = load_benchmark_results(rawFile);
nSolvers = sum([cfg.solvers.enabled]);
expectedCases = 15 * 2 * 3 - 2;     % guess8 excluded for P02 at both n
assert(height(T) == expectedCases * nSolvers, ...
    'Pipeline produced %d rows, expected %d.', height(T), expectedCases * nSolvers);
assert(~any(T.status == "implementation_error"), ...
    'Implementation error in pipeline: %s', strjoin(T.message(T.status == "implementation_error"), ' | '));
% resume must not duplicate rows
benchmark_controller(cfg);
T2 = load_benchmark_results(rawFile);
assert(height(T2) == height(T), 'Resume duplicated runs.');
export_all_outputs(cfg);
required = {fullfile('tables', 'summary_by_solver.csv'), ...
    fullfile('tables', 'summary_by_solver.tex'), fullfile('tables', 'detail_general.tex'), ...
    fullfile('tables', 'detail_spd.tex'), fullfile('tables', 'results_wide.csv'), ...
    fullfile('tables', 'restart_statistics.csv'), fullfile('tables', 'failure_breakdown.csv'), ...
    fullfile('profiles', 'profile_general_iterations.pdf'), ...
    fullfile('profiles', 'profile_spd_elapsed_seconds.csv'), ...
    fullfile('profiles', 'profile_spd_capped_function_evaluations.png'), ...
    fullfile('histories_plots', 'history_P13_n10_guess1.pdf'), 'protocol.txt'};
for i = 1:numel(required)
    assert(exist(fullfile(cfg.outputDirectory, required{i}), 'file') == 2, ...
        'Missing output %s.', required{i});
end
end
