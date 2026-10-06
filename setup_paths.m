function root = setup_paths()
%SETUP_PATHS Add the framework folders to the MATLAB path.
%   ROOT = SETUP_PATHS() adds solvers/, problems/, utils/ and analysis/
%   (relative to this file) to the path and returns the framework root.
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'solvers'), fullfile(root, 'problems'), ...
    fullfile(root, 'utils'), fullfile(root, 'analysis'));
end
