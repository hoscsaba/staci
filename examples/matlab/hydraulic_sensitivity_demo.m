function sensitivity = hydraulic_sensitivity_demo()
%HYDRAULIC_SENSITIVITY_DEMO Pressure and flow sensitivity to one diameter.
exampleDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(fileparts(exampleDirectory));
addpath(fullfile(repositoryDirectory, 'matlab'));

network = StaciModel(fullfile(repositoryDirectory, 'tests', ...
    'anytown_1med.spr'));
pipeId = network.linkIds("pipe");
pipeId = pipeId(1);
status = network.solveHydraulics();
assert(status.converged);

sensitivity = network.hydraulicSensitivity(pipeId, "diameter_m");
pressureTable = table(string(sensitivity.node_ids(:)), ...
    sensitivity.node_pressure_head_derivative(:), ...
    'VariableNames', {'node_id', 'd_pressure_head_d_diameter'});
fprintf('Sensitivity parameter: %s diameter [m]\n', pipeId);
disp(pressureTable);
end
