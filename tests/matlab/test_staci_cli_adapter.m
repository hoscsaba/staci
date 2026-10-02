function test_staci_cli_adapter(binary)
%TEST_STACI_CLI_ADAPTER Exercise the external-process example, without MEX.
% addpath('tests/matlab'); test_staci_cli_adapter('/absolute/path/build/staci')
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    oldPath = path;
    pathCleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
    addpath(fullfile(root, 'examples', 'integration'));
    work = [tempname ' MATLAB jobs'];
    mkdir(work);
    filesCleanup = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>
    network = fullfile(root, 'tests', 'epanet_eps_smoke.inp');
    for epsMode = [false true]
        job = run_staci_cli(binary, network, work, 'EPS', epsMode);
        assert(strcmp(job.state, 'success') && job.exit_code == 0);
        assert(~isempty(job.run_id) && isfile(job.result_file));
        result = jsondecode(fileread(job.result_file));
        assert(isstruct(result));
    end
    badInput = fullfile(work, 'invalid input.inp');
    fid = fopen(badInput, 'w');
    assert(fid >= 0);
    fprintf(fid, '[JUNCTIONS]\nJ bad-number\n[END]\n');
    fclose(fid);
    job = run_staci_cli(binary, badInput, work);
    assert(strcmp(job.state, 'failed') && job.exit_code == 2);
    assert(isempty(job.result_file) && ~isempty(job.diagnostics));
    job = run_staci_cli(binary, network, work, 'EPS', true, 'Timeout', 1e-9);
    assert(strcmp(job.state, 'timeout') && isempty(job.exit_code));
    assert(isempty(job.result_file));
    fprintf('PASS: MATLAB CLI steady/EPS, invalid input, timeout and paths with spaces.\n');
end
