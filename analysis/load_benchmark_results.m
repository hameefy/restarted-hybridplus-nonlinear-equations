function T = load_benchmark_results(rawFile)
%LOAD_BENCHMARK_RESULTS Read raw_runs.csv into a typed table.
%
%   T = LOAD_BENCHMARK_RESULTS(RAWFILE) reads the file written by
%   BENCHMARK_CONTROLLER. Text columns are strings; SUCCESS is logical.
%   If a (case_id, solver) pair occurs more than once (e.g. after a manual
%   rerun), only the LAST row is kept and a warning is issued.

if ~exist(rawFile, 'file')
    error('bm:io', 'Raw results file not found: %s', rawFile);
end
opts = detectImportOptions(rawFile, 'Delimiter', ',', 'TextType', 'string');
textColumns = {'case_id', 'group', 'problem', 'start', 'solver', 'status', ...
    'timestamp', 'message'};
opts = setvartype(opts, intersect(textColumns, opts.VariableNames), 'string');
T = readtable(rawFile, opts);
if height(T) == 0
    error('bm:io', 'Raw results file %s contains no runs.', rawFile);
end
T.success = logical(T.success);

key = T.case_id + "|" + T.solver;
[~, lastIdx] = unique(key, 'last');
if numel(lastIdx) < height(T)
    warning('bm:duplicateRuns', ...
        '%d duplicate (case, solver) rows found; keeping the last occurrence.', ...
        height(T) - numel(lastIdx));
    T = T(sort(lastIdx), :);
end
end
