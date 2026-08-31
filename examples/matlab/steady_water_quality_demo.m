function result = steady_water_quality_demo(showFigure)
%STEADY_WATER_QUALITY_DEMO Solve hydraulics and asymptotic water age.
if nargin < 1
    showFigure = true;
end
exampleDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(fileparts(exampleDirectory));
addpath(fullfile(repositoryDirectory, 'matlab'));

network = StaciModel(fullfile(repositoryDirectory, 'tests', ...
    'anytown_1med.spr'));
status = network.solveHydraulics();
assert(status.converged, 'STACI:QualityDemo', ...
    'Hydraulics must converge before the water-quality calculation.');
quality = network.solveSteadyWaterAge();
nodes = network.nodeResults();

result = table(string(quality.node_ids(:)), nodes.demand_m3s(:), ...
    quality.node_water_age_s(:), quality.node_water_age_s(:) / 3600, ...
    'VariableNames', {'node_id', 'demand_m3s', 'water_age_s', ...
                      'water_age_h'});
disp(result);

if showFigure
    figure('Color', 'white', 'Name', 'STACI steady water age');
    bar(result.water_age_h);
    ylabel('Steady water age [h]');
    xlabel('Node index');
    grid on;
    title('Anytown asymptotic water age');
end
end
