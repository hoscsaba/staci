function comparison = modify_network_in_memory()
%MODIFY_NETWORK_IN_MEMORY Change pipe diameters and re-solve without files.
exampleDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(fileparts(exampleDirectory));
addpath(fullfile(repositoryDirectory, 'matlab'));

network = StaciModel(fullfile(repositoryDirectory, 'tests', ...
    'anytown_1med.spr'));
pipeIds = network.linkIds("pipe");
initialDiameters = network.getLinkProperty(pipeIds, "diameter_m");

initialStatus = network.solveHydraulics();
assert(initialStatus.converged);
initialNodes = network.nodeResults();

% This is a vectorized MEX call; the complete model remains in C++ memory.
network.setLinkProperty(pipeIds, "diameter_m", 1.10 * initialDiameters);
modifiedStatus = network.solveHydraulics();
assert(modifiedStatus.converged);
modifiedNodes = network.nodeResults();

comparison = table(string(initialNodes.id(:)), ...
    initialNodes.pressure_head_m(:), modifiedNodes.pressure_head_m(:), ...
    modifiedNodes.pressure_head_m(:) - initialNodes.pressure_head_m(:), ...
    'VariableNames', {'node_id', 'initial_pressure_head_m', ...
                      'modified_pressure_head_m', 'change_m'});
disp(comparison);

% Restore the original diameters in memory. No SPR/INP file was modified.
network.setLinkProperty(pipeIds, "diameter_m", initialDiameters);
end
