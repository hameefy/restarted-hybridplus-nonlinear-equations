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
%   Plot options (cfg.analysis):
%     plotTitles  true: title "<metric>: <subset> (<cases>)" above the plot;
%                 false (default): no title, for figures with a caption.
%     logTauAxis  true (default): logarithmic tau axis with ticks at powers
%                 of two, labelled in tau; false: linear axis (or log2(tau)
%                 if log2Scale is true).

T = load_benchmark_results(fullfile(resultsDir, 'raw_runs.csv'));
outDir = fullfile(resultsDir, 'profiles');
if ~exist(outDir, 'dir'), mkdir(outDir); end
solverNames = ordered_solvers(T, cfg);
subsets = benchmark_subsets(T, cfg);
styles = {'-', '--', '-.', ':', '-', '--', '-.'};
markers = {'none', 'none', 'none', 'none', 'o', 's', 'd'};
showTitles = get_option(cfg.analysis, 'plotTitles', false);
logTau = get_option(cfg.analysis, 'logTauAxis', true);

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
        if ~logTau && cfg.analysis.log2Scale
            xValues = log2(tauPlot);
            xLabel = 'log_2(\tau)';
        end

        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 640 440]);
        ax = axes(fig); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on'); %#ok<LAXES> one new figure per profile
        set(ax, 'FontSize', 13);
        for s = 1:numel(solverNames)
            [xs, ys] = stairs(xValues, rhoPlot(:, s));   % staircase coordinates
            plot(ax, xs, ys, 'LineStyle', styles{mod(s-1, 7)+1}, ...
                'Marker', markers{mod(s-1, 7)+1}, 'MarkerIndices', ...
                unique(round(linspace(1, numel(xs), 8))), 'LineWidth', 2);
        end
        if logTau
            set_log_tau_axis(ax, tauMax);
        end
        xlim(ax, [xValues(1), xValues(end)]);
        ylim(ax, [0, 1.02]);
        xlabel(ax, xLabel); ylabel(ax, '\rho_s(\tau)');
        if showTitles
            title(ax, sprintf('%s: %s (%d cases)', label, subsets(b).description, ...
                numel(prof.caseIds)), 'Interpreter', 'none', 'FontWeight', 'normal');
        end
        legend(ax, solverNames, 'Location', 'best', 'Interpreter', 'none', ...
            'FontSize', 12);
        base = fullfile(outDir, sprintf('profile_%s_%s', subsets(b).name, metric));
        exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector', 'Padding', 6);
        exportgraphics(fig, [base '.png'], 'Resolution', 200, 'Padding', 6);
        close(fig);

        curve = array2table([prof.tau, prof.rho], 'VariableNames', ...
            [{'tau'}, solverNames]);
        writetable(curve, [base '.csv']);
    end
end
end

function value = get_option(s, name, default)
value = default;
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
end
end

function set_log_tau_axis(ax, tauMax)
% Logarithmic tau axis with at most about nine ticks at powers of two,
% labelled by their tau values (1, 2, 4, ...), or as 2^k beyond 1024.
set(ax, 'XScale', 'log');
kMax = max(1, ceil(log2(tauMax)));
step = max(1, ceil((kMax + 1) / 9));
k = 0:step:kMax;
labels = cell(size(k));
for i = 1:numel(k)
    if k(i) <= 10
        labels{i} = sprintf('%d', 2^k(i));
    else
        labels{i} = sprintf('2^{%d}', k(i));
    end
end
set(ax, 'XTick', 2.^k, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex', ...
    'XMinorTick', 'off', 'XMinorGrid', 'off');
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
