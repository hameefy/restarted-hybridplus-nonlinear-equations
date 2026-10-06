function export_all_outputs(cfg)
%EXPORT_ALL_OUTPUTS Regenerate every table and figure from raw_runs.csv.
%
%   EXPORT_ALL_OUTPUTS(CFG) produces the tables (GENERATE_TABLES), the
%   performance profiles (GENERATE_PERFORMANCE_PROFILES) and the
%   convergence-history plots for cfg.analysis.historyCases
%   (PLOT_CONVERGENCE_HISTORY) in cfg.outputDirectory. It never reruns a
%   solver, so it can be called any number of times after the benchmark.

resultsDir = cfg.outputDirectory;
generate_tables(resultsDir, cfg);
generate_performance_profiles(resultsDir, cfg);
histDir = fullfile(resultsDir, 'histories');
figDir = fullfile(resultsDir, 'histories_plots');
if ~exist(figDir, 'dir'), mkdir(figDir); end
for i = 1:numel(cfg.analysis.historyCases)
    caseId = cfg.analysis.historyCases{i};
    if exist(fullfile(histDir, [caseId '.mat']), 'file')
        plot_convergence_history(resultsDir, caseId, {}, ...
            fullfile(figDir, ['history_' caseId '.pdf']));
    else
        warning('bm:analysis', 'No history saved for %s; plot skipped.', caseId);
    end
end
fprintf('Tables, profiles and history plots written to %s\n', resultsDir);
end
