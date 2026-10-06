function profile = compute_performance_profile(T, metric, solverNames, floorValue)
%COMPUTE_PERFORMANCE_PROFILE Dolan-More performance profile from raw runs.
%
%   PROFILE = COMPUTE_PERFORMANCE_PROFILE(T, METRIC, SOLVERNAMES, FLOOR)
%   uses the rows of table T (see LOAD_BENCHMARK_RESULTS) for the solvers
%   in SOLVERNAMES and the cost column METRIC ('iterations',
%   'function_evaluations' or 'elapsed_seconds').
%
%   Definition (Dolan and More, Math. Program. 91 (2002) 201-213):
%     t_{p,s} = cost of solver s on case p if the run succeeded, else Inf;
%     r_{p,s} = t_{p,s} / min_s t_{p,s};
%     rho_s(tau) = |{p : r_{p,s} <= tau}| / |P|.
%   Conventions (stated in the protocol):
%   * The case set P is every case_id in T; its size is the denominator.
%     A case that no solver solves stays in P with r = Inf for all.
%   * A (case, solver) pair missing from T counts as a failure (Inf).
%   * Costs are floored at FLOOR (e.g. 1 iteration, 1e-6 s) so that zero
%     costs (a start that already satisfies the tolerance, or a timer
%     resolution of zero) do not produce 0/0 ratios.
%
%   PROFILE fields: solvers, caseIds, cost (|P|-by-S), ratio (|P|-by-S),
%   tau (column of breakpoints, starting at 1), rho (numel(tau)-by-S),
%   solvedFraction (1-by-S, rho at tau = Inf).

if ~ismember(metric, T.Properties.VariableNames)
    error('bm:analysis', 'Unknown metric column "%s".', metric);
end
solverNames = cellstr(solverNames);
caseIds = unique(T.case_id);
nc = numel(caseIds);
ns = numel(solverNames);
if nc == 0
    error('bm:analysis', 'No cases available for the performance profile.');
end
cost = Inf(nc, ns);
[~, caseIndex] = ismember(T.case_id, caseIds);
for s = 1:ns
    rows = find(T.solver == solverNames{s});
    for r = rows'
        if T.success(r) && isfinite(T.(metric)(r))
            cost(caseIndex(r), s) = max(T.(metric)(r), floorValue);
        end
    end
end
best = min(cost, [], 2);
ratio = cost ./ best;
ratio(~isfinite(best), :) = Inf;          % nobody solved the case

finiteRatios = ratio(isfinite(ratio));
tau = unique([1; finiteRatios(:)]);
rho = zeros(numel(tau), ns);
for s = 1:ns
    rs = sort(ratio(:, s));
    for k = 1:numel(tau)
        rho(k, s) = sum(rs <= tau(k)) / nc;
    end
end
profile = struct('metric', metric, 'solvers', {solverNames}, ...
    'caseIds', {caseIds}, 'cost', cost, 'ratio', ratio, 'tau', tau, ...
    'rho', rho, 'solvedFraction', sum(isfinite(ratio), 1) / nc);
end
