function tests = test_staci_matlab_interface
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(fileparts(testDirectory));
addpath(fullfile(repositoryDirectory, 'matlab'));
addpath(fullfile(repositoryDirectory, 'examples', 'matlab'));
if exist('staci_mex', 'file') ~= 3
    build_staci_mex();
end
testCase.TestData.repositoryDirectory = repositoryDirectory;
testCase.TestData.networkFile = fullfile(repositoryDirectory, 'tests', ...
    'anytown_1med.spr');
end

function testLifecycleAndIntrospection(testCase)
model = StaciModel(testCase.TestData.networkFile);
verifyEqual(testCase, numel(model.nodeIds()), 23);
verifyEqual(testCase, numel(model.linkIds("pipe")), 41);
verifyEqual(testCase, string(staci_mex('version')), ...
    "STACI MATLAB API 1.0");
links = model.linkIds();
verifyError(testCase, @() model.getLinkProperty( ...
    links(1), "not_a_property"), 'STACI:MatlabAPI');
model.release();
verifyError(testCase, @() model.nodeIds(), 'STACI:InvalidHandle');
end

function testVectorPropertyRoundTripAndHydraulics(testCase)
model = StaciModel(testCase.TestData.networkFile);
pipes = model.linkIds("pipe");
initial = model.getLinkProperty(pipes, "diameter_m");
modified = 1.05 * initial;
model.setLinkProperty(pipes, "diameter_m", modified);
verifyEqual(testCase, model.getLinkProperty(pipes, "diameter_m"), ...
    modified, 'AbsTol', 1e-12);
status = model.solveHydraulics();
verifyTrue(testCase, status.converged);
nodes = model.nodeResults();
links = model.linkResults();
verifySize(testCase, nodes.pressure_head_m, [23, 1]);
verifySize(testCase, links.flow_rate_m3s, [44, 1]);
verifyTrue(testCase, all(isfinite(nodes.pressure_head_m)));
end

function testInpLoadAndHydraulics(testCase)
networkFile = fullfile(testCase.TestData.repositoryDirectory, 'tests', ...
    'epanet_smoke.inp');
model = StaciModel(networkFile);
verifyGreaterThan(testCase, numel(model.nodeIds()), 0);
verifyGreaterThan(testCase, numel(model.linkIds()), 0);
status = model.solveHydraulics();
verifyTrue(testCase, status.converged);
verifyTrue(testCase, all(isfinite(model.nodeResults().pressure_head_m)));
end

function testSteadyWaterAge(testCase)
model = StaciModel(testCase.TestData.networkFile);
verifyTrue(testCase, model.solveHydraulics().converged);
quality = model.solveSteadyWaterAge();
verifySize(testCase, quality.node_water_age_s, [23, 1]);
verifyTrue(testCase, all(isfinite(quality.node_water_age_s)));
verifyGreaterThanOrEqual(testCase, quality.node_water_age_s, zeros(23, 1));
end

function testDiameterSensitivityAgainstCentralDifference(testCase)
model = StaciModel(testCase.TestData.networkFile);
pipe = model.linkIds("pipe");
pipe = pipe(1);
diameter = model.getLinkProperty(pipe, "diameter_m");
verifyTrue(testCase, model.solveHydraulics().converged);
sensitivity = model.hydraulicSensitivity(pipe, "diameter_m");

delta = 1e-3 * diameter;
model.setLinkProperty(pipe, "diameter_m", diameter + delta);
model.resetHydraulicState();
verifyTrue(testCase, model.solveHydraulics().converged);
plus = model.nodeResults();
model.setLinkProperty(pipe, "diameter_m", diameter - delta);
model.resetHydraulicState();
verifyTrue(testCase, model.solveHydraulics().converged);
minus = model.nodeResults();
finiteDifference = (plus.pressure_head_m - minus.pressure_head_m) / (2 * delta);
relativeError = norm(sensitivity.node_pressure_head_derivative - ...
    finiteDifference) / norm(finiteDifference);
verifyLessThan(testCase, relativeError, 0.05);
end

function testExamples(testCase)
hydraulics = getting_started_hydraulics(false);
verifyTrue(testCase, hydraulics.status.converged);
comparison = modify_network_in_memory();
verifyGreaterThan(testCase, max(abs(comparison.change_m)), 1e-8);
quality = steady_water_quality_demo(false);
verifyEqual(testCase, height(quality), 23);
sensitivity = hydraulic_sensitivity_demo();
verifySize(testCase, sensitivity.node_pressure_head_derivative, [23, 1]);
end
