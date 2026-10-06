function bm_validate_config(cfg)
%BM_VALIDATE_CONFIG Check a benchmark configuration before any run starts.
%
%   BM_VALIDATE_CONFIG(CFG) raises an error with identifier 'bm:config' if
%   CFG (see BENCHMARK_CONFIG) is incomplete or inconsistent, so that a
%   long experiment never fails after hours because of a setup mistake.

required = {'dimensions', 'problemIds', 'startIds', 'seed', 'protocol', ...
    'solvers', 'timingRepeats', 'warmUp', 'resume', 'stopOnError', ...
    'saveHistories', 'verbose', 'outputDirectory', 'analysis'};
missing = setdiff(required, fieldnames(cfg));
if ~isempty(missing)
    error('bm:config', 'Configuration is missing field(s): %s.', strjoin(missing, ', '));
end
if ~isfield(cfg, 'initialStepRule') || ~ischar(cfg.initialStepRule) || ...
        ~ismember(cfg.initialStepRule, {'unit', 'secant', 'fd'})
    error('bm:config', ['cfg.initialStepRule must be ''unit'', ''secant'' or ', ...
        '''fd'' (see benchmark_config.m).']);
end
if isempty(cfg.dimensions) || any(cfg.dimensions < 1) || ...
        any(cfg.dimensions ~= floor(cfg.dimensions))
    error('bm:config', 'cfg.dimensions must contain positive integers.');
end
registry = problem_registry();
unknownProblems = setdiff(cfg.problemIds, [registry.id]);
if ~isempty(unknownProblems)
    error('bm:config', 'Unknown problem id(s): %s.', mat2str(unknownProblems));
end
if isempty(cfg.startIds) || any(~ismember(cfg.startIds, 1:10))
    error('bm:config', 'cfg.startIds must be a subset of 1:10.');
end
if ~(isscalar(cfg.seed) && cfg.seed >= 0 && cfg.seed == floor(cfg.seed))
    error('bm:config', 'cfg.seed must be a nonnegative integer.');
end
if ~(isscalar(cfg.timingRepeats) && cfg.timingRepeats >= 1 && ...
        cfg.timingRepeats == floor(cfg.timingRepeats))
    error('bm:config', 'cfg.timingRepeats must be a positive integer.');
end
protocolFields = fieldnames(bm_common_defaults());
if ~isempty(setxor(fieldnames(cfg.protocol), protocolFields))
    error('bm:config', 'cfg.protocol must have exactly the fields: %s.', ...
        strjoin(protocolFields', ', '));
end
if isempty(cfg.solvers) || ~any([cfg.solvers.enabled])
    error('bm:config', 'At least one solver must be enabled.');
end
names = {cfg.solvers.name};
if numel(unique(names)) ~= numel(names)
    error('bm:config', 'Solver names must be unique.');
end
for s = 1:numel(cfg.solvers)
    solver = cfg.solvers(s);
    if ~isa(solver.fun, 'function_handle')
        error('bm:config', 'Solver %s: fun must be a function handle.', solver.name);
    end
    if ~isvarname(solver.name)
        error('bm:config', 'Solver name "%s" must be a valid MATLAB identifier.', ...
            solver.name);
    end
    overlap = intersect(fieldnames(solver.options), protocolFields);
    if ~isempty(overlap)
        error('bm:config', ['Solver %s overrides protocol option(s) %s; ', ...
            'protocol options must be identical for all solvers.'], ...
            solver.name, strjoin(overlap', ', '));
    end
end
end
