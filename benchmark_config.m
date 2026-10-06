function cfg = benchmark_config()
%BENCHMARK_CONFIG All experimental settings in one visible place.
%
%   CFG = BENCHMARK_CONFIG() returns the configuration of the benchmark:
%   restarted hybridplus, hybridplus_cr, dprpplus and gnbfgs on P1-P15,
%   n in {6500, 11500, 550000, 1000000}, guess1-guess10. Edit this file (or
%   modify the returned structure before calling BENCHMARK_CONTROLLER) to
%   change the experiment. Nothing in the solvers or the controller
%   hard-codes any of these choices.
%
%   Parameters:
%   * Stopping rule: success iff ||F(x_k)|| <= 1e-6 within 1000 iterations
%     and 2000 residual evaluations; otherwise "F".
%   * Starts: guess1-guess10; guess8 excluded for P2.
%   * restarted_hybridplus / hybridplus_cr: ell = 0.95, u = 1 (satisfy
%     the compatibility condition for P11-P15), t = 0.95, sigma = 1e-4,
%     rho = 0.5, tau_k = 2^{-k}; secant initial trial clipped to
%     [1e-6, 1e6].
%   * dprpplus: theta = 1, rho = 0.2 (paper); sigma = 0.01,
%     tau_k = exp(-(k+1)^(k+1)), ell = 1e-10, u = 1e10 (not stated in the
%     paper; chosen by the authors).
%   * gnbfgs: r = 0.1, rho = sqrt(0.9), sigma1 = sigma2 = 1e-5,
%     lambda_{-1} = 0.01, eps_k = 1/(k+1)^2, B_0 = I (Li-Fukushima 1999);
%     dense matrix, so runs are capped at n <= 11500.

% ------------------------------------------------------------------
% Hybridplus parameters of this benchmark
% ------------------------------------------------------------------
cfg.initialStepRule = 'secant';
cfg.t = 0.95;

% ------------------------------------------------------------------
% Benchmark suite
% ------------------------------------------------------------------
cfg.dimensions = [6500, 11500, 550000, 1000000];
cfg.problemIds = 1:15;            % see PROBLEM_REGISTRY
cfg.startIds = 1:10;              % guess1..guess10 (see PROBLEM_STARTS)
cfg.seed = 20260930;              % random starts use RandStream(seed + n)

% ------------------------------------------------------------------
% Protocol (identical for every solver)
% ------------------------------------------------------------------
cfg.protocol = struct( ...
    'tolerance', 1e-6, ...
    'maxIterations', 1000, ...
    'maxFunctionEvaluations', 2000, ...
    'recordHistory', true);

% ------------------------------------------------------------------
% Solvers: name, function handle, enabled flag, solver-specific options.
% Every solver has the signature  result = solver(problem, x0, options).
% ------------------------------------------------------------------
cfg.solvers = struct('name', {}, 'fun', {}, 'enabled', {}, 'options', {});
cfg.solvers(end+1) = struct('name', 'restarted_hybridplus', ...
    'fun', @restarted_hybridplus, 'enabled', true, 'options', struct( ...
    'ell', 0.95, 'u', 1, 't', cfg.t, 'sigma', 1e-4, 'rho', 0.5, ...
    'tauFcn', @(k) 2.^(-k), 'initialStep', cfg.initialStepRule, ...
    'stepMin', 1e-6, 'stepMax', 1e6));
cfg.solvers(end+1) = struct('name', 'hybridplus_cr', ...
    'fun', @hybridplus_cr, 'enabled', true, 'options', struct( ...
    't', cfg.t, 'sigma', 1e-4, 'rho', 0.5, 'tauFcn', @(k) 2.^(-k), ...
    'initialStep', cfg.initialStepRule, 'stepMin', 1e-6, 'stepMax', 1e6));
cfg.solvers(end+1) = struct('name', 'dprpplus', ...
    'fun', @dprpplus, 'enabled', true, 'options', struct( ...
    'ell', 1e-10, 'u', 1e10, 'theta', 1, 'rho', 0.2, 'sigma', 0.01, ...
    'tauFcn', @(k) exp(-(k + 1).^(k + 1))));
cfg.solvers(end+1) = struct('name', 'gnbfgs', ...
    'fun', @gnbfgs, 'enabled', true, 'options', struct( ...
    'r', 0.1, 'rho', sqrt(0.9), 'sigma1', 1e-5, 'sigma2', 1e-5, ...
    'lambdaMinus1', 0.01, 'epsFcn', @(k) 1 ./ (k + 1).^2, ...
    'maxDimension', 11500));

% ------------------------------------------------------------------
% Execution
% ------------------------------------------------------------------
cfg.timingRepeats = 1;      % >1: rerun each case, report the median time
cfg.warmUp = true;          % one untimed call per solver before the suite
cfg.resume = true;          % skip (case, solver) pairs already in raw_runs.csv
cfg.stopOnError = false;    % false: log implementation errors and continue
cfg.saveHistories = true;   % per-case residual histories (for plots)
cfg.verbose = true;
cfg.outputDirectory = fullfile(fileparts(mfilename('fullpath')), 'results');

% ------------------------------------------------------------------
% Analysis (performance profiles and tables)
% ------------------------------------------------------------------
cfg.analysis.metrics = {'iterations', 'function_evaluations', 'elapsed_seconds'};
cfg.analysis.metricLabels = {'Iter', 'FuncEvals', 'ElapsedSec'};
cfg.analysis.metricFloor = struct('iterations', 1, ...
    'function_evaluations', 1, 'elapsed_seconds', 1e-6);
cfg.analysis.tauMax = [];          % [] = automatic (largest finite ratio)
cfg.analysis.log2Scale = false;    % true: plot against log2(tau)
cfg.analysis.referenceSolver = 'restarted_hybridplus';
cfg.analysis.historyCases = {'P11_n6500_guess1', 'P05_n6500_guess1'};
end
