%MAIN_RUN_BENCHMARKS Main driver: full benchmark, then all tables/figures.
%
%   Run from any folder:
%       run('<path>/restarted-hybridplus-nonlinear-equations/main_run_benchmarks.m')
%   or, with the repository folder as the current folder:
%       main_run_benchmarks
%
%   Steps
%   1. Add the framework to the path.
%   2. Load and validate the configuration (benchmark_config.m).
%   3. Record the environment and every setting (results/protocol.txt).
%   4. Run all enabled solvers on all cases (results/raw_runs.csv,
%      results/histories/*.mat). Interrupted runs resume automatically.
%   5. Build tables, performance profiles and history plots from the raw
%      file only (results/tables, results/profiles, results/histories_plots).
%
%   Order of work:
%     a. run_self_checks          software tests (about a minute)
%     b. main_run_benchmarks      restarted_hybridplus, hybridplus_cr,
%                                 dprpplus, gnbfgs -> results/
%   To run a subset, edit the configuration after loading it, e.g.
%       cfg.dimensions = [6500 11500];
%       cfg.problemIds = 11:15;

setup_paths();
cfg = benchmark_config();
bm_validate_config(cfg);

fprintf('Output directory: %s\n', cfg.outputDirectory);
fprintf('Enabled solvers: %s\n', strjoin({cfg.solvers([cfg.solvers.enabled]).name}, ', '));
write_protocol(cfg);
rawFile = benchmark_controller(cfg);
fprintf('Raw results: %s\n', rawFile);
export_all_outputs(cfg);
