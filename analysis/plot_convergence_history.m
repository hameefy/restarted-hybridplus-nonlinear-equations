function plot_convergence_history(resultsDir, caseId, solverNames, outFile)
%PLOT_CONVERGENCE_HISTORY Residual norm versus iteration for one case.
%
%   PLOT_CONVERGENCE_HISTORY(RESULTSDIR, CASEID) plots ||F(x_k)|| against k
%   (logarithmic vertical axis) for every solver stored in
%   RESULTSDIR/histories/<CASEID>.mat, e.g. CASEID = 'P11_n1500_guess1'.
%   PLOT_CONVERGENCE_HISTORY(..., SOLVERNAMES) restricts to the given cell
%   array of solver names. PLOT_CONVERGENCE_HISTORY(..., OUTFILE) saves the
%   figure (vector PDF if OUTFILE ends in .pdf, otherwise raster) and
%   closes it; without OUTFILE the figure stays open.
%   Restart iterations of the hybridplus variants are marked with circles.
%   The tolerance line is taken from raw_runs.csv when available.

if nargin < 3 || isempty(solverNames), solverNames = {}; end
if nargin < 4, outFile = ''; end
file = fullfile(resultsDir, 'histories', [caseId '.mat']);
if ~exist(file, 'file')
    error('bm:io', 'No history file for case %s (%s).', caseId, file);
end
data = load(file, 'caseHistory');
available = fieldnames(data.caseHistory);
if isempty(solverNames)
    solverNames = available;
end
solverNames = cellstr(solverNames);
missing = setdiff(solverNames, available);
if ~isempty(missing)
    error('bm:analysis', 'Case %s has no history for: %s.', caseId, strjoin(missing, ', '));
end

visible = 'on';
if ~isempty(outFile), visible = 'off'; end
fig = figure('Visible', visible, 'Color', 'w', 'Position', [100 100 640 440]);
ax = axes(fig); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
set(ax, 'YScale', 'log');
styles = {'-', '--', '-.', ':'};
handles = gobjects(0);
for s = 1:numel(solverNames)
    h = data.caseHistory.(solverNames{s});
    if ~h.enabled
        warning('bm:analysis', 'History for %s was not recorded; skipped.', solverNames{s});
        continue;
    end
    k = (0:numel(h.residualNorm) - 1)';
    res = max(h.residualNorm, realmin);
    handles(end+1) = plot(ax, k, res, 'LineStyle', styles{mod(s-1, 4)+1}, ...
        'LineWidth', 1.6, 'DisplayName', solverNames{s}); %#ok<AGROW>
    if any(h.restarted)
        plot(ax, k(h.restarted), res(h.restarted), 'o', 'Color', ...
            handles(end).Color, 'HandleVisibility', 'off');
    end
end
rawFile = fullfile(resultsDir, 'raw_runs.csv');
if exist(rawFile, 'file')
    T = load_benchmark_results(rawFile);
    tol = T.tolerance(find(T.case_id == caseId, 1));
    if ~isempty(tol)
        yline(ax, tol, ':k', 'tolerance', 'HandleVisibility', 'off');
    end
end
xlabel(ax, 'Iteration k'); ylabel(ax, '||F(x_k)||');
title(ax, strrep(caseId, '_', ' '), 'FontWeight', 'normal');
legend(ax, handles, 'Location', 'northeast', 'Interpreter', 'none');
if ~isempty(outFile)
    [~, ~, ext] = fileparts(outFile);
    if strcmpi(ext, '.pdf')
        exportgraphics(fig, outFile, 'ContentType', 'vector');
    else
        exportgraphics(fig, outFile, 'Resolution', 200);
    end
    close(fig);
end
end
