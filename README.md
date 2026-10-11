# Restarted hybridplus for nonlinear equations — MATLAB code and benchmark results

This repository contains the MATLAB code and the complete benchmark results for the paper

> H. Mohammad and H. Choi, *A restarted Hestenes–Stiefel/Polak–Ribière–Polyak hybrid method with bounded diagonal scaling for nonlinear systems of equations* (submitted).

It compares four solvers on the same benchmark suite: restarted hybridplus, hybridplus_cr, dprpplus and gnbfgs. Every table and figure of the paper's numerical section is generated from a single raw log, `results/raw_runs.csv`.

**Tested with:** MATLAB R2026a on Windows 11. Checks performed:
- `checkcode` reports 0 messages in all 43 files;
- `run_self_checks` passes all 9 checks.

The framework uses `exportgraphics`, `groupcounts`, `yline` and string tables, so it should run on R2020a or later. Only R2026a has actually been tested.

**Licence:** MIT (see `LICENSE`).

## 1. Quick start

**The results reported in the paper are in `results/`.** To regenerate their tables and figures without rerunning anything:

```matlab
cd <path>/restarted-hybridplus-nonlinear-equations
setup_paths; cfg = benchmark_config(); export_all_outputs(cfg)
```

**To reproduce the experiments from scratch**, first rename or move the folder `results/`. Otherwise `main_run_benchmarks` resumes from the existing `results/raw_runs.csv` and skips every run that is already recorded there.

```matlab
cd <path>/restarted-hybridplus-nonlinear-equations
run_self_checks          % about a minute; must print "All self checks passed."
main_run_benchmarks      % all four solvers -> results/
```

If the suite is interrupted, run `main_run_benchmarks` again. Runs already in `raw_runs.csv` are skipped (`cfg.resume = true`).

**Expected run time:** several hours on a desktop machine. Most of it is gnbfgs, whose dense n×n matrix makes it by far the slowest solver at n = 6500 and n = 11500.

To run a subset (for example to check timing first):

```matlab
setup_paths; cfg = benchmark_config();
cfg.dimensions = 6500; cfg.problemIds = 11:15;
cfg.outputDirectory = fullfile(pwd, 'results_pilot');
write_protocol(cfg); benchmark_controller(cfg); export_all_outputs(cfg);
```

## 2. File tree

```
restarted-hybridplus-nonlinear-equations/
  LICENSE, README.md
  results/                  the results reported in the paper (see §6)
  main_run_benchmarks.m     main driver script
  benchmark_config.m        ALL settings (suite, protocol, solver parameters, analysis)
  benchmark_controller.m    paired run loop, timing, raw CSV logging, resume
  run_self_checks.m         software tests (no benchmark data)
  setup_paths.m
  solvers/  restarted_hybridplus.m  hybridplus_cr.m  hybridplus_core.m (shared by the two)
            dprpplus.m  gnbfgs.m
  problems/ problem_P01.m ... problem_P15.m  problem_registry.m  problem_starts.m
  utils/    ResidualOracle.m (counted, budget-enforcing F)  bm_*.m helpers
  analysis/ load_benchmark_results.m  compute_performance_profile.m
            generate_performance_profiles.m  generate_tables.m
            plot_convergence_history.m  benchmark_subsets.m
            write_protocol.m  export_all_outputs.m
```

## 3. Uniform interfaces

**Solver:** `result = solver(problem, x0, options)`
- `problem.F` is a function handle `x -> F(x)` taking and returning column vectors.
- `options` contains the protocol fields (`tolerance`, `maxIterations`, `maxFunctionEvaluations`, `recordHistory`) plus the solver's own parameters. **Unknown fields are rejected**, so a misspelled parameter cannot be silently ignored.
- `result` has the fields `solver, x, finalResidual, status, success, message, iterations, functionEvaluations, invalidEvaluations, jacobianEvaluations (=0), restarts, history`.

**Residual evaluations.** Every solver evaluates F only through `ResidualOracle`. Every call counts: the initial point, every rejected positive or negative trial, and the two auxiliary gnbfgs evaluations. The oracle refuses the call after the budget is exhausted, so a run can never exceed 2000 evaluations. A value that is not finite or not real (for example log of a negative number in P02) is reported as invalid, and the trial is rejected.

**Success rule:** a run succeeds iff ‖F(x_k)‖₂ ≤ 10⁻⁶ within 1000 iterations and 2000 residual evaluations. Every other run is a failure, reported as "F" in the tables.

## 4. Adding a problem or a solver

- **Problem:** write `problems/problem_P16.m` with `F = problem_P16(x)`. Add one `add(...)` line in `problem_registry.m`, giving its group, minimal dimension, excluded starts, LaTeX definition and notes. Then add `16` to `cfg.problemIds`.
- **Solver:** write `solvers/mysolver.m` with the interface above. Use `bm_merge_options`, `bm_validate_common` (which builds the oracle), `bm_history` and `bm_make_result`; `dprpplus.m` is the shortest template. Then append an entry to `cfg.solvers` in `benchmark_config.m`.

## 5. Status taxonomy

| status | meaning | counts as |
|---|---|---|
| `success` | ‖F‖ ≤ tolerance | success |
| `max_iterations` | 1000 accepted iterations reached | failure |
| `max_function_evaluations` | 2000 residual evaluations reached | failure |
| `line_search_failure` | step length underflowed to 0 | failure |
| `numerical_error` | NaN or Inf in the search direction | failure |
| `invalid_initial_point` | F(x₀) not finite or not real | failure |
| `invalid_trial_point` | gnbfgs: g(x_k + λ_{k−1}g_k) not finite, so p_k is undefined | failure |
| `not_run_dimension_cap` | gnbfgs not attempted because n > 11500 (dense n×n matrix) | failure (see §7) |
| `implementation_error` | exception in a solver or problem (message logged) | failure |

## 6. Outputs (in `results/`)

| file | content |
|---|---|
| `raw_runs.csv` | one row per (case, solver): status, success, iterations, function_evaluations, invalid_evaluations, restarts, final_residual, elapsed_seconds (tic/toc), cpu_seconds (cputime), budgets, timestamp, message |
| `histories/<case>.mat` | ‖F(x_k)‖, step lengths and restart flags per iteration for every solver |
| `protocol.txt`, `config_snapshot.mat` | MATLAB version, OS, CPU, threads, and every parameter |
| `tables/summary_by_solver.{csv,tex}` | solved counts, failures by type, wins, medians on commonly solved cases, restarts |
| `tables/detail_general.tex`, `tables/detail_spd.tex` | per-case longtables (Iter, FuncEvals, ElapsedSec, FinalResidual; "F" = failure, "NR" = not run because of the gnbfgs dimension cap) |
| `tables/results_wide.csv`, `restart_statistics.csv`, `failure_breakdown.csv` | supplementary data |
| `profiles/profile_<subset>_<metric>.{pdf,png,csv}` | Dolan–Moré profiles; subsets: general, spd, all, general_capped, spd_capped |
| `histories_plots/history_<case>.pdf` | convergence histories for `cfg.analysis.historyCases` |

## 7. Parameters that must be reported

All of these are set in `benchmark_config.m` and recorded in `protocol.txt`.

**restarted_hybridplus / hybridplus_cr:** ℓ = 0.95, u = 1, t = 0.95, σ = 10⁻⁴, ρ = 0.5, τ_k = 2⁻ᵏ.
- The initial trial step is the secant rule λ⁰ = |Fᵀd| / (μ‖d‖²), with μ = |sᵀy| / ‖s‖² and λ⁰₀ = 1. It is clipped to [λ_min, λ_max] = [10⁻⁶, 10⁶] and needs no extra evaluation.
- (ℓ, u) satisfy the compatibility condition for P11–P15.
- hybridplus_cr is restarted hybridplus with B_k = I.
- The code also offers the initial-trial options `unit` (λ⁰ = 1) and `fd` (finite-difference, one extra evaluation per iteration). The benchmark does not use them.

**dprpplus:** θ = 1 and ρ = 0.2 are from the paper. σ = 0.01, τ_k = exp(−(k+1)^(k+1)), ℓ = 10⁻¹⁰ and u = 10¹⁰ were chosen by the authors, because the paper does not state them.
- In double precision τ_k = 0 for k ≥ 4.
- The coefficient's denominator is implemented **as printed**, max{‖F_{k−1}‖², ‖d_{k−1}‖}, using the norm of d_{k−1}.

**gnbfgs:** r = 0.1, ρ = √0.9, σ₁ = σ₂ = 10⁻⁵, λ₋₁ = 0.01, B₀ = I, from the paper.
- The paper prints ε_k = k⁻², which is undefined at k = 0. The framework uses ε_k = (k+1)⁻².
- The inverse BFGS form is used; in exact arithmetic it is equivalent to the paper's form.
- **Dense n×n storage:** runs are attempted only for n ≤ 11500 (about 1 GB, with a peak of about 3 GB during the update). At n = 550000 the matrix would need about 2.4 TB. On the uncapped subsets these cases count as failures; use the `*_capped` profiles for a like-for-like comparison.

**Problems and starts** (definitions in `problems/` and in the manuscript):
- P1: F₁ = e^{x₁} − 1 and Fᵢ = e^{xᵢ} + xᵢ − 1 for i = 2, …, n.
- P6: Fₙ = xₙ(xₙ₋₁² + xₙ²), as in Li and Fukushima (1999).
- P10 uses h = 2.
- P2 excludes guess8 (domain xᵢ > −1).
- guess6ᵢ = (n + 1 − 2i)/(n + 1), i = 1, …, n. All ten starts are distinct.
- **P9** computes Fᵢ = xᵢ − 4/(c·n·(δᵀx)·Σⱼ δⱼ⁻¹), with δⱼ = (j − ½)/n and c = 0.9: "a function derived from the Chandrasekhar H-equation".

**Random starts:** guess9 and guess10 are drawn by `RandStream('mt19937ar','Seed', 20260930 + n)` (guess9 first, then guess10, which is an independent draw with its sign flipped). They are identical for all problems and solvers.

## 8. Run counts and resources

- Dimensions: n ∈ {6500, 11500, 550000, 1000000}.
- Per solver: P1–P10 give 10 × 4 × 10 − 4 = **396** cases; P11–P15 give **200**; the total is **596**. With 4 solvers that is 2384 runs.
- gnbfgs is attempted on 298 cases (n ≤ 11500); the other 298 are `not_run_dimension_cap`.
- All problems cost O(n) per evaluation. At n = 10⁶ a residual vector is 8 MB.
- Timing runs single-threaded or multi-threaded according to MATLAB defaults, and the thread count is recorded. For stable timings, close other applications.
