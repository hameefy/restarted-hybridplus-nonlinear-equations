function write_protocol(cfg)
%WRITE_PROTOCOL Record environment and every setting of the run.
%
%   WRITE_PROTOCOL(CFG) writes <cfg.outputDirectory>/protocol.txt and
%   <cfg.outputDirectory>/config_snapshot.mat. Call it immediately before
%   BENCHMARK_CONTROLLER so that the recorded environment is the one used.
%   If protocol.txt already exists (resumed run), a timestamped copy is
%   written instead of overwriting it.

outDir = cfg.outputDirectory;
if ~exist(outDir, 'dir'), mkdir(outDir); end
file = fullfile(outDir, 'protocol.txt');
if exist(file, 'file')
    file = fullfile(outDir, sprintf('protocol_%s.txt', ...
        char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
end
fid = fopen(file, 'w');
if fid < 0, error('bm:io', 'Cannot write %s.', file); end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Restarted hybridplus benchmark protocol\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'MATLAB: %s\n', version);
fprintf(fid, 'Computer: %s\n', computer);
fprintf(fid, 'Operating system: %s\n', system_description());
fprintf(fid, 'CPU cores (physical, MATLAB): %d\n', feature('numcores'));
fprintf(fid, 'maxNumCompThreads: %d\n', maxNumCompThreads);
fprintf(fid, 'CPU: %s\n', cpu_description());
fprintf(fid, '\n[Suite]\n');
fprintf(fid, 'Dimensions: %s\n', mat2str(cfg.dimensions));
fprintf(fid, 'Problems: %s\n', mat2str(cfg.problemIds));
fprintf(fid, 'Starts: %s (guess8 excluded for P02)\n', mat2str(cfg.startIds));
fprintf(fid, 'Random starts: RandStream(''mt19937ar'', ''Seed'', %d + n); guess9 then guess10\n', cfg.seed);
fprintf(fid, '\n[Protocol, identical for all solvers]\n');
fprintf(fid, 'Success iff ||F(x_k)||_2 <= %g\n', cfg.protocol.tolerance);
fprintf(fid, 'Max iterations: %d; max residual evaluations: %d\n', ...
    cfg.protocol.maxIterations, cfg.protocol.maxFunctionEvaluations);
fprintf(fid, 'Every call of F counts (initial point, rejected trials, probes).\n');
fprintf(fid, 'Timing: tic/toc wall time and cputime of the solver call; repeats = %d (median); warm-up = %d\n', ...
    cfg.timingRepeats, cfg.warmUp);
fprintf(fid, '\n[Solvers]\n');
for s = 1:numel(cfg.solvers)
    fprintf(fid, '%s (enabled = %d): %s\n', cfg.solvers(s).name, ...
        cfg.solvers(s).enabled, func2str(cfg.solvers(s).fun));
    names = fieldnames(cfg.solvers(s).options);
    for i = 1:numel(names)
        v = cfg.solvers(s).options.(names{i});
        if isa(v, 'function_handle')
            v = func2str(v);
        else
            v = num2str(v, 16);
        end
        fprintf(fid, '    %s = %s\n', names{i}, v);
    end
end
fprintf(fid, '\n[Analysis]\n');
fprintf(fid, 'Profile metrics: %s\n', strjoin(cfg.analysis.metrics, ', '));
floors = cfg.analysis.metricFloor;
names = fieldnames(floors);
for i = 1:numel(names)
    fprintf(fid, 'Cost floor %s = %g\n', names{i}, floors.(names{i}));
end
fprintf(fid, 'Failures and missing runs have ratio Inf; every case stays in the denominator.\n');

save(fullfile(outDir, 'config_snapshot.mat'), 'cfg');
end

function s = system_description()
if ispc
    [~, s] = system('ver');
elseif ismac
    [~, s] = system('sw_vers -productVersion');
else
    [~, s] = system('uname -srm');
end
s = strtrim(s);
end

function s = cpu_description()
s = 'unknown';
try
    if ispc
        [status, out] = system('wmic cpu get name');
        if status == 0
            lines = strtrim(splitlines(string(out)));
            lines = lines(lines ~= "" & lines ~= "Name");
            if ~isempty(lines), s = char(lines(1)); end
        end
    elseif ismac
        [~, out] = system('sysctl -n machdep.cpu.brand_string');
        s = strtrim(out);
    else
        [~, out] = system('grep -m1 "model name" /proc/cpuinfo');
        s = strtrim(out);
    end
catch
    s = 'unknown';
end
end
