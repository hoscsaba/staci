function results = run_matlab_tests()
%RUN_MATLAB_TESTS Run the optional MATLAB/MEX developer tests.
testDirectory = fileparts(mfilename('fullpath'));
results = runtests(testDirectory, 'IncludeSubfolders', true);
disp(results);
assert(all([results.Passed]), 'STACI:MatlabTestsFailed', ...
    'One or more STACI MATLAB tests failed.');
end
