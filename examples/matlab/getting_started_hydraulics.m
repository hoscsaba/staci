function result = getting_started_hydraulics(showFigure)
%GETTING_STARTED_HYDRAULICS Load, solve and inspect a STACI network in memory.
if nargin < 1
    showFigure = true;
end
exampleDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(fileparts(exampleDirectory));
addpath(fullfile(repositoryDirectory, 'matlab'));

networkFile = fullfile(repositoryDirectory, 'tests', 'anytown_1med.spr');
network = StaciModel(networkFile);
fprintf('Loaded %d nodes and %d links.\n', ...
    numel(network.nodeIds()), numel(network.linkIds()));

status = network.solveHydraulics();
assert(status.converged, 'STACI:GettingStarted', ...
    'The hydraulic calculation did not converge.');
nodes = network.nodeTable();
links = network.linkTable();

fprintf('Pressure-head range: %.3f ... %.3f m\n', ...
    min(nodes.pressure_head_m), max(nodes.pressure_head_m));
fprintf('Maximum absolute link flow: %.6f m3/s\n', ...
    max(abs(links.flow_rate_m3s)));

if showFigure
    figure('Color', 'white', 'Name', 'STACI MATLAB hydraulic results');
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    nexttile;
    bar(nodes.pressure_head_m);
    ylabel('Pressure head [m]');
    xlabel('Node index');
    grid on;
    title('Anytown node pressures');
    nexttile;
    bar(links.flow_rate_m3s);
    ylabel('Flow rate [m^3/s]');
    xlabel('Link index');
    grid on;
    title('Anytown signed link flows');
end

result = struct('status', status, 'nodes', nodes, 'links', links);
end
