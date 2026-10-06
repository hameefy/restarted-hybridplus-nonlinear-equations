function rawFile = benchmark_controller(cfg)
%BENCHMARK_CONTROLLER Run every (problem, dimension, start, solver) case.
%
%   RAWFILE = BENCHMARK_CONTROLLER(CFG) runs the benchmark defined by CFG
%   (see BENCHMARK_CONFIG) and appends one row per run to
%   <cfg.outputDirectory>/raw_runs.csv, which is returned in RAWFILE. The
%   raw file is the single source for every table and figure.
%
%   Pairing: for each problem and dimension, every enabled solver receives
%   exactly the same initial points (PROBLEM_STARTS with cfg.seed).
%   Timing: wall-clock time of the solver call (tic/toc) and process CPU
%   time (cputime), excluding problem construction and logging. With
%   cfg.timingRepeats > 1 the call is repeated and the median times are
%   reported; counters and status are taken from the first call (all
%   solvers are deterministic).
%   Failures: an exception thrown by a solver or a problem function is
%   recorded with status 'implementation_error' (and the run continues
%   unless cfg.stopOnError is true). No failed run is ever dropped.
%   Resume: with cfg.resume = true, (case, solver) pairs already present in
%   raw_runs.csv are skipped, so an interrupted suite can be restarted.

bm_validate_config(cfg);
outDir = cfg.outputDirectory;
if ~exist(outDir, 'dir'), mkdir(outDir); end
histDir = fullfile(outDir, 'histories');
if cfg.saveHistories && ~exist(histDir, 'dir'), mkdir(histDir); end

rawFile = fullfile(outDir, 'raw_runs.csv');
columns = {'run_id', 'case_id', 'group', 'problem', 'dimension', 'start', ...
    'solver', 'status', 'success', 'iterations', 'function_evaluations', ...
    'invalid_evaluations', 'restarts', 'final_residual', 'elapsed_seconds', ...
    'cpu_seconds', 'tolerance', 'max_iterations', 'max_function_evaluations', ...
    'timestamp', 'message'};
done = containers.Map('KeyType', 'char', 'ValueType', 'logical');
runId = 0;
if exist(rawFile, 'file') && cfg.resume
    previous = readtable(rawFile, 'TextType', 'string', 'Delimiter', ',');
    for r = 1:height(previous)
        done(char(previous.case_id(r) + "|" + previous.solver(r))) = true;
    end
    if height(previous) > 0, runId = max(previous.run_id); end
    fid = fopen(rawFile, 'a');
elseif exist(rawFile, 'file') && ~cfg.resume
    error('bm:outputExists', ['%s already exists and cfg.resume is false. ', ...
        'Move or delete it, or set cfg.resume = true.'], rawFile);
else
    fid = fopen(rawFile, 'w');
    fprintf(fid, '%s\n', strjoin(columns, ','));
end
if fid < 0
    error('bm:io', 'Cannot open %s for writing.', rawFile);
end
cleanup = onCleanup(@() fclose(fid));

registry = problem_registry();
solvers = cfg.solvers([cfg.solvers.enabled]);
if cfg.warmUp
    warm_up(solvers, cfg);
end

nCases = 0;
for problemId = cfg.problemIds
    entry = registry([registry.id] == problemId);
    for n = cfg.dimensions
        if n < entry.minDimension, continue; end
        starts = problem_starts(n, cfg.seed);
        for startId = cfg.startIds
            if ismember(startId, entry.excludedStarts), continue; end
            start = starts(startId);
            caseId = sprintf('%s_n%d_%s', entry.name, n, start.name);
            problem = struct('id', entry.id, 'name', entry.name, ...
                'group', entry.group, 'dimension', n, 'F', entry.F);
            caseHistory = struct();
            for s = 1:numel(solvers)
                solver = solvers(s);
                key = [caseId '|' solver.name];
                if isKey(done, key), continue; end
                runId = runId + 1;
                options = merge_structs(cfg.protocol, solver.options);
                [result, wall, cpu] = timed_call(solver, problem, start.x0, ...
                    options, cfg);
                row = {runId, caseId, entry.group, entry.name, n, start.name, ...
                    solver.name, result.status, double(result.success), ...
                    result.iterations, result.functionEvaluations, ...
                    result.invalidEvaluations, result.restarts, ...
                    result.finalResidual, wall, cpu, cfg.protocol.tolerance, ...
                    cfg.protocol.maxIterations, cfg.protocol.maxFunctionEvaluations, ...
                    char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')), ...
                    result.message};
                write_row(fid, row);
                if cfg.saveHistories
                    caseHistory.(solver.name) = result.history;
                end
                if cfg.verbose
                    fprintf('%6d  %-24s %-22s %-34s it=%4d fe=%5d res=%9.2e t=%8.3fs\n', ...
                        runId, caseId, solver.name, result.status, result.iterations, ...
                        result.functionEvaluations, result.finalResidual, wall);
                end
                clear result
            end
            if cfg.saveHistories && ~isempty(fieldnames(caseHistory))
                save_history(fullfile(histDir, [caseId '.mat']), caseHistory);
            end
            nCases = nCases + 1;
        end
    end
end
if cfg.verbose
    fprintf('Benchmark finished: %d problem-dimension-start cases processed.\n', nCases);
end
end

% ======================================================================
function [result, wall, cpu] = timed_call(solver, problem, x0, options, cfg)
walls = zeros(cfg.timingRepeats, 1);
cpus = zeros(cfg.timingRepeats, 1);
result = [];
for rep = 1:cfg.timingRepeats
    try
        c0 = cputime;
        t0 = tic;
        current = solver.fun(problem, x0, options);
        walls(rep) = toc(t0);
        cpus(rep) = cputime - c0;
    catch err
        if cfg.stopOnError
            rethrow(err);
        end
        result = failure_result(solver.name, err);
        wall = NaN; cpu = NaN;
        return;
    end
    if rep == 1
        result = current;
    end
    clear current
end
check_result(result, solver.name);
result.x = [];                       % do not keep large final iterates
wall = median(walls);
cpu = median(cpus);
end

function result = failure_result(solverName, err)
msg = sprintf('%s: %s', err.identifier, err.message);
result = struct('solver', solverName, 'x', [], 'finalResidual', Inf, ...
    'status', 'implementation_error', 'success', false, 'message', msg, ...
    'iterations', 0, 'functionEvaluations', NaN, 'invalidEvaluations', NaN, ...
    'jacobianEvaluations', 0, 'restarts', 0, ...
    'history', struct('enabled', false, 'residualNorm', NaN, ...
    'stepLength', NaN, 'restarted', false));
end

function check_result(result, solverName)
required = {'status', 'success', 'iterations', 'functionEvaluations', ...
    'invalidEvaluations', 'restarts', 'finalResidual', 'message', 'history'};
missing = setdiff(required, fieldnames(result));
if ~isempty(missing)
    error('bm:resultSchema', '%s returned a result without field(s): %s.', ...
        solverName, strjoin(missing, ', '));
end
end

function out = merge_structs(a, b)
out = a;
names = fieldnames(b);
for i = 1:numel(names)
    out.(names{i}) = b.(names{i});
end
end

function write_row(fid, row)
parts = cell(size(row));
for i = 1:numel(row)
    v = row{i};
    if ischar(v) || isstring(v)
        v = char(v);
        v = strrep(v, '"', '''');
        v = regexprep(v, '[\r\n]+', ' ');
        parts{i} = ['"' v '"'];
    elseif islogical(v)
        parts{i} = sprintf('%d', v);
    elseif isnumeric(v) && isscalar(v) && v == floor(v) && abs(v) < 1e15
        parts{i} = sprintf('%d', v);
    else
        parts{i} = sprintf('%.16g', v);
    end
end
fprintf(fid, '%s\n', strjoin(parts, ','));
end

function save_history(file, caseHistory)
if exist(file, 'file')
    old = load(file, 'caseHistory');
    names = fieldnames(caseHistory);
    for i = 1:numel(names)
        old.caseHistory.(names{i}) = caseHistory.(names{i});
    end
    caseHistory = old.caseHistory;
end
save(file, 'caseHistory');
end

function warm_up(solvers, cfg)
% One untimed call per solver on a small SPD problem, to exclude
% just-in-time compilation and first-call overheads from the timings.
problem = struct('id', 13, 'name', 'P13', 'group', 'spd', 'dimension', 10, ...
    'F', @problem_P13);
x0 = ones(10, 1);
for s = 1:numel(solvers)
    options = merge_structs(cfg.protocol, solvers(s).options);
    solvers(s).fun(problem, x0, options);
end
end
