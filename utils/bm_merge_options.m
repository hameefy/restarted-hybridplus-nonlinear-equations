function options = bm_merge_options(options, defaults, solverName)
%BM_MERGE_OPTIONS Fill missing option fields from DEFAULTS; reject unknown ones.
%
%   OPTIONS = BM_MERGE_OPTIONS(OPTIONS, DEFAULTS, SOLVERNAME) returns a
%   structure containing every field of DEFAULTS. Fields present in OPTIONS
%   override the defaults. A field of OPTIONS that does not appear in
%   DEFAULTS raises an error, so that a misspelled parameter can never be
%   silently ignored.

if nargin < 1 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error('bm:invalidOptions', '%s: options must be a scalar structure.', solverName);
end
given = fieldnames(options);
known = fieldnames(defaults);
unknown = setdiff(given, known);
if ~isempty(unknown)
    error('bm:invalidOptions', '%s: unknown option(s): %s.', solverName, ...
        strjoin(unknown', ', '));
end
for i = 1:numel(known)
    name = known{i};
    if ~isfield(options, name) || isempty(options.(name))
        options.(name) = defaults.(name);
    end
end
end
