function bm_require_scalar(value, name, solverName, lower, upper, closedLower, closedUpper)
%BM_REQUIRE_SCALAR Validate a real scalar parameter against an interval.
%
%   BM_REQUIRE_SCALAR(VALUE, NAME, SOLVERNAME, LOWER, UPPER) requires
%   LOWER < VALUE < UPPER. CLOSEDLOWER (default false) and CLOSEDUPPER
%   (default false) make the corresponding end point admissible.

if nargin < 6 || isempty(closedLower), closedLower = false; end
if nargin < 7 || isempty(closedUpper), closedUpper = false; end
ok = isnumeric(value) && isscalar(value) && isreal(value) && ~isnan(value);
if ok
    if closedLower, ok = value >= lower; else, ok = value > lower; end
end
if ok
    if closedUpper, ok = value <= upper; else, ok = value < upper; end
end
if ~ok
    lb = '('; if closedLower, lb = '['; end
    rb = ')'; if closedUpper, rb = ']'; end
    error('bm:invalidOptions', '%s: %s must be a real scalar in %s%g, %g%s.', ...
        solverName, name, lb, lower, upper, rb);
end
end
