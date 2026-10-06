function generate_performance_profiles(resultsDir, cfg)
%GENERATE_PERFORMANCE_PROFILES Dolan-More profiles for every subset/metric.
%
%   GENERATE_PERFORMANCE_PROFILES(RESULTSDIR, CFG) reads
%   RESULTSDIR/raw_runs.csv and writes, for each subset of BENCHMARK_SUBSETS
%   and each metric in cfg.analysis.metrics,
%     RESULTSDIR/profiles/profile_<subset>_<metric>.pdf / .png
%     RESULTSDIR/profiles/profile_<subset>_<metric>.csv   (tau, rho per solver)
%   The solvers are those present in the raw file, in the order of
%   cfg.solvers. Figures are created invisibly and closed after saving.

T = load_benchmark_results(fullfile(resultsDir, 'raw_runs.csv'));
outDir = fullfile(resultsDir, 'profiles');
if ~exist(outDir, 'dir'), mkdir(outDir); end
solverNames = ordered_solvers(T, cfg);
subsets = benchmark_subsets(T, cfg);
styles = {'-', '--', '-.', ':', '-', '--', '-.'};
markers = {'none', 'none', 'none', 'none', 'o', 's', 'd'};

for b = 1:numel(subsets)
    Ts = T(subsets(b).mask, :);
    for m = 1:numel(cfg.analysis.metrics)
        metric = cfg.analysis.metrics{m};
        label = cfg.analysis.metricLabels{m};
        prof = compute_performance_profile(Ts, metric, solverNames, ...
            cfg.analysis.metricFloor.(metric));
        tauMax = cfg.analysis.tauMax;
        if isempty(tauMax)
            tauMax = max([prof.tau; 2]) * 1.05;
        end
        tauPlot = [prof.tau; tauMax];
        rhoPlot = [prof.rho; prof.rho(end, :)];
        xValues = tauPlot;
        xLabel = '\tau';
        if cfg.analysis.log2Scale
            xValues = log2(tauPlot);
            xLabel = 'log_2(\tau)';
        end

        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 640 440]);
        ax = axes(fig); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on'); %#ok<LAXES> one new figure per profile
        for s = 1:numel(solverNames)
            [xs, ys] = stairs(xValues, rhoPlot(:, s));   % staircase coordinates
            plot(ax, xs, ys, 'LineStyle', styles{mod(s-1, 7)+1}, ...
                'Marker', markers{mod(s-1, 7)+1}, 'MarkerIndices', ...
                unique(round(linspace(1, numel(xs), 8))), 'LineWidth', 1.6);
        end
        xlim(ax, [xValues(1), xValues(end)]);
        ylim(ax, [0, 1.02]);
        xlabel(ax, xLabel); ylabel(ax, '\rho_s(\tau)');
        title(ax, sprintf('%s: %s (%d cases)', label, subsets(b).description, ...
            numel(prof.caseIds)), 'Interpreter', 'none', 'FontWeight', 'normal');
        legend(ax, solverNames, 'Location', 'southeast', 'Interpreter', 'none');
        base = fullfile(outDir, sprintf('profile_%s_%s', subsets(b).name, metric));
        exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
        exportgraphics(fig, [base '.png'], 'Resolution', 200);
        close(fig);

        curve = array2table([prof.tau, prof.rho], 'VariableNames', ...
            [{'tau'}, solverNames]);
        writetable(curve, [base '.csv']);
    end
end
end

function names = ordered_solvers(T, cfg)
present = unique(T.solver);
names = {};
for s = 1:numel(cfg.solvers)
    if any(present == cfg.solvers(s).name)
        names{end+1} = cfg.solvers(s).name; %#ok<AGROW>
    end
end
extra = setdiff(cellstr(present), names);
names = [names, extra(:)'];
end
