function mexFile = build_staci_mex(buildDirectory)
%BUILD_STACI_MEX Build the general in-memory STACI MATLAB interface.

scriptDirectory = fileparts(mfilename('fullpath'));
sourceDirectory = fileparts(scriptDirectory);
if nargin < 1 || isempty(buildDirectory) || strlength(string(buildDirectory)) == 0
    buildDirectory = fullfile(sourceDirectory, 'build-matlab');
end

if exist('staci_mex', 'file') == 3
    clear staci_mex
end

cmakeExecutable = findCMake();
configureCommand = sprintf([ ...
    '"%s" -S "%s" -B "%s" ' ...
    '-DSTACI_BUILD_MATLAB_MEX=ON ' ...
    '-DSTACI_BUILD_OPTIMIZERS=OFF ' ...
    '-DSTACI_ENABLE_HDF5=OFF ' ...
    '-DBUILD_TESTING=OFF ' ...
    '-DMatlab_ROOT_DIR="%s" ' ...
    '-DCMAKE_BUILD_TYPE=Release'], ...
    cmakeExecutable, sourceDirectory, buildDirectory, matlabroot);
runCommand(configureCommand, 'CMake configuration');

buildCommand = sprintf( ...
    '"%s" --build "%s" --config Release --target staci_mex --parallel', ...
    cmakeExecutable, buildDirectory);
runCommand(buildCommand, 'MEX compilation');

mexDirectories = {fullfile(buildDirectory, 'matlab'), ...
                  fullfile(buildDirectory, 'Matlab')};
for index = 1:numel(mexDirectories)
    if isfolder(mexDirectories{index})
        addpath(mexDirectories{index});
    end
end
mexCandidates = dir(fullfile(buildDirectory, '**', ...
    ['staci_mex.', mexext]));
for index = 1:numel(mexCandidates)
    addpath(mexCandidates(index).folder);
end
mexFile = which('staci_mex');
if isempty(mexFile)
    error('STACI:MatlabBuild', ...
        'The build completed but staci_mex is not on the MATLAB path.');
end
fprintf('STACI MATLAB interface: %s\n', mexFile);
end

function executable = findCMake()
configured = strtrim(getenv('CMAKE_COMMAND'));
candidates = {configured, ...
    '/opt/homebrew/bin/cmake', ...
    '/usr/local/bin/cmake', ...
    '/Applications/CMake.app/Contents/bin/cmake', ...
    fullfile(getenv('ProgramFiles'), 'CMake', 'bin', 'cmake.exe')};
for index = 1:numel(candidates)
    candidate = candidates{index};
    if ~isempty(candidate) && isfile(candidate)
        executable = candidate;
        return
    end
end
[status, output] = system(ternary(ispc, 'where cmake', 'command -v cmake'));
if status == 0
    discovered = regexp(strtrim(output), '\r?\n', 'split');
    for index = 1:numel(discovered)
        if isfile(strtrim(discovered{index}))
            executable = strtrim(discovered{index});
            return
        end
    end
end
error('STACI:MissingCMake', [ ...
    'CMake was not found. Install it (for example: brew install cmake), ' ...
    'or set CMAKE_COMMAND to the full cmake executable path.']);
end

function value = ternary(condition, trueValue, falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end

function runCommand(command, description)
fprintf('%s...\n%s\n', description, command);
[status, output] = system(command);
fprintf('%s', output);
if status ~= 0
    error('STACI:MatlabBuild', '%s failed with exit status %d.', ...
        description, status);
end
end
