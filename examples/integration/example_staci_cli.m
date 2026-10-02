function jobs = example_staci_cli(binary, network, jobRoot)
%EXAMPLE_STACI_CLI Steady hydraulics and EPS calls with result/error handling.
% example_staci_cli('/absolute/path/staci', '/absolute/path/network.inp', 'jobs')
% Use an INP with a configured duration for the EPS example.

    jobs.steady = run_staci_cli(binary, network, jobRoot);
    displayDiagnostics(jobs.steady);
    if strcmp(jobs.steady.state, 'success') && ~isempty(jobs.steady.result_file)
        hydraulics = jsondecode(fileread(jobs.steady.result_file));
        disp(hydraulics);
    end

    jobs.eps = run_staci_cli(binary, network, jobRoot, 'EPS', true, 'Timeout', 300);
    displayDiagnostics(jobs.eps);
    if strcmp(jobs.eps.state, 'success') && ~isempty(jobs.eps.result_file)
        metadata = jsondecode(fileread(jobs.eps.result_file));
        nodes = readtable(fullfile(jobs.eps.job_dir, 'result-nodes.csv'));
        disp(metadata);
        disp(nodes(1:min(5, height(nodes)), :));
    elseif strcmp(jobs.eps.state, 'partial')
        fprintf('Partial EPS output: inspect frame convergence before use.\n');
    end
end

function displayDiagnostics(job)
    fprintf('STACI state: %s; job directory: %s\n', job.state, job.job_dir);
    for k = 1:numel(job.diagnostics)
        r = job.diagnostics{k};
        if isfield(r, 'severity') && ismember(r.severity, {'warning', 'error'})
            fprintf(2, '%s: %s\n', r.code, r.message);
        end
    end
end
