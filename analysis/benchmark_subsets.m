function subsets = benchmark_subsets(T, cfg)
%BENCHMARK_SUBSETS Case subsets used for profiles and summary tables.
%
%   SUBSETS = BENCHMARK_SUBSETS(T, CFG) returns a struct array with fields
%   name, description and mask (logical rows of T):
%     general        P1-P10, all dimensions
%     spd            P11-P15, all dimensions
%     all            every case
%     general_capped P1-P10 with n <= dimension cap of a capped solver
%     spd_capped     P11-P15 with n <= that cap
%   The *_capped subsets exist only if an enabled solver has an option
%   maxDimension < max(cfg.dimensions) (gnbfgs by default). On those
%   subsets every solver was actually attempted on every case, so they give
%   the like-for-like comparison including that solver. On the uncapped
%   subsets, 'not_run_dimension_cap' counts as a failure; this must be
%   stated whenever those profiles are reported.

subsets = struct('name', {}, 'description', {}, 'mask', {});
subsets(end+1) = struct('name', 'general', ...
    'description', 'General problems P1-P10', 'mask', T.group == "general");
subsets(end+1) = struct('name', 'spd', ...
    'description', 'SPD problems P11-P15', 'mask', T.group == "spd");
subsets(end+1) = struct('name', 'all', ...
    'description', 'All problems P1-P15', 'mask', true(height(T), 1));

cap = Inf;
for s = find([cfg.solvers.enabled])
    if isfield(cfg.solvers(s).options, 'maxDimension')
        cap = min(cap, cfg.solvers(s).options.maxDimension);
    end
end
if isfinite(cap) && cap < max(cfg.dimensions)
    subsets(end+1) = struct('name', 'general_capped', 'description', ...
        sprintf('General problems P1-P10, n <= %g', cap), ...
        'mask', T.group == "general" & T.dimension <= cap);
    subsets(end+1) = struct('name', 'spd_capped', 'description', ...
        sprintf('SPD problems P11-P15, n <= %g', cap), ...
        'mask', T.group == "spd" & T.dimension <= cap);
end
subsets = subsets(arrayfun(@(s) any(s.mask), subsets));
end
