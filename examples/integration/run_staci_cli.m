function job = run_staci_cli(binary, network, jobRoot, varargin)
%RUN_STACI_CLI Run the STACI executable in an isolated job directory.
% job = run_staci_cli(binary, network, jobRoot, 'EPS', true, 'Timeout', 300)
% Requires MATLAB with JVM support; no MEX module or Python is needed.
% Diagnostics is a cell array because JSONL records have different fields.
% States: success, partial, failed, timeout. exit_code is [] after timeout.
% Process startup/filesystem errors raise MATLAB exceptions.

    p = inputParser;
    addParameter(p, 'EPS', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'Timeout', 120, @(x) isnumeric(x) && isscalar(x) && ...
        isfinite(x) && x > 0);
    parse(p, varargin{:});
    assert(usejava('jvm'), 'STACI:Adapter:JVM', 'MATLAB JVM support is required.');
    binary = absolutePath(binary);
    network = absolutePath(network);
    jobRoot = absolutePath(jobRoot);
    assert(isfile(binary), 'STACI:Adapter:Binary', 'Executable not found: %s', binary);
    assert(isfile(network), 'STACI:Adapter:Input', 'Input not found: %s', network);
    ensureDirectory(jobRoot);
    work = tempname(jobRoot);
    ensureDirectory(work);
    [~, ~, ext] = fileparts(network);
    copied = fullfile(work, ['network' lower(ext)]);
    [ok, message] = copyfile(network, copied);
    assert(ok, 'STACI:Adapter:Copy', '%s', message);
    diagnosticsPath = fullfile(work, 'diagnostics.jsonl');
    args = {binary, '--diagnostics-file', diagnosticsPath};
    if p.Results.EPS
        args = [args, {'--epanet-eps', copied, '-o', fullfile(work, 'result')}];
        result = fullfile(work, 'result.meta.json');
    else
        args = [args, {'-s', copied}];
        result = [copied '.hydraulics.json'];
    end

    % Pass an argument list directly: spaces and shell characters are literal.
    command = javaObject('java.util.ArrayList');
    for k = 1:numel(args)
        command.add(javaObject('java.lang.String', args{k}));
    end
    builder = javaObject('java.lang.ProcessBuilder', command);
    builder.directory(javaObject('java.io.File', work));
    builder.redirectErrorStream(true);
    builder.redirectOutput(javaObject('java.io.File', fullfile(work, 'console.log')));
    process = builder.start();
    processCleanup = onCleanup(@() stopProcess(process)); %#ok<NASGU>
    started = tic;
    timedOut = false;
    while process.isAlive()
        if toc(started) >= p.Results.Timeout
            timedOut = true;
            stopProcess(process);
            break;
        end
        pause(min(0.05, p.Results.Timeout));
    end
    exitCode = [];
    if ~timedOut
        exitCode = double(process.exitValue());
    end

    records = {};
    if isfile(diagnosticsPath)
        lines = regexp(fileread(diagnosticsPath), '\r?\n', 'split');
        for k = 1:numel(lines)
            if isempty(strtrim(lines{k})), continue; end
            try
                record = jsondecode(lines{k});
                if ~isstruct(record) || ~isscalar(record)
                    error('STACI:Adapter:Record', 'Expected a JSON object.');
                end
            catch
                record = struct('severity', 'warning', ...
                    'code', 'ADAPTER.INCOMPLETE_LOG', ...
                    'message', 'Incomplete or invalid diagnostic record.');
            end
            records{end + 1} = record; %#ok<AGROW>
        end
    end
    runId = '';
    for k = 1:numel(records)
        r = records{k};
        if isfield(r, 'event') && strcmp(r.event, 'run_start') && isfield(r, 'run_id')
            runId = r.run_id;
        end
    end
    complete = false;
    for k = numel(records):-1:1
        r = records{k};
        if ~isempty(runId) && isfield(r, 'event') && strcmp(r.event, 'run_end') && ...
                isfield(r, 'run_id') && strcmp(r.run_id, runId)
            complete = isfield(r, 'exit_code') && isequal(r.exit_code, exitCode);
            break;
        end
    end
    state = 'failed';
    if timedOut
        state = 'timeout';
    elseif complete && exitCode == 0
        state = 'success';
    elseif complete && exitCode == 3
        state = 'partial';
    end
    if ~ismember(state, {'success', 'partial'}) || ~isfile(result)
        result = '';
    end
    job = struct('state', state, 'exit_code', exitCode, 'run_id', runId, ...
        'job_dir', work, 'diagnostics', {records}, 'result_file', result);
end

function path = absolutePath(path)
    assert((ischar(path) && isrow(path)) || (isstring(path) && isscalar(path)), ...
        'STACI:Adapter:Path', 'Paths must be character vectors or scalar strings.');
    file = javaObject('java.io.File', char(path));
    path = char(file.getCanonicalPath());
end

function ensureDirectory(path)
    [ok, message] = mkdir(path);
    assert(ok, 'STACI:Adapter:Directory', '%s', message);
end

function stopProcess(process)
    if process.isAlive()
        process.destroyForcibly();
        process.waitFor();
    end
end
