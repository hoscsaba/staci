function result = optimize_anytown_residence_time(networkFile, userOptions)
%OPTIMIZE_ANYTOWN_RESIDENCE_TIME Continuous pipe-diameter GA demonstration.
%
% result = optimize_anytown_residence_time()
% result = optimize_anytown_residence_time(networkFile)
% result = optimize_anytown_residence_time(networkFile, gaOptions)
%
% Objective: demand-weighted mean steady water age [s].
% Variables: every Cso diameter, continuously bounded to 0.01...1.00 m.
% Constraints: every calculated nodal pressure head is in 0...60 m.
%
% StaciModel keeps the parsed STACI model in memory. No SPR/INP or result
% file is written during GA objective evaluations.

scriptDirectory = fileparts(mfilename('fullpath'));
sourceDirectory = fileparts(fileparts(scriptDirectory));
if nargin < 1 || isempty(networkFile) || strlength(string(networkFile)) == 0
    networkFile = fullfile(sourceDirectory, 'tests', 'anytown_1med.spr');
end
addpath(fullfile(sourceDirectory, 'matlab'));
model = StaciModel(networkFile);
pipeIds = model.linkIds("pipe");
nodeIds = model.nodeIds();
originalDiameters = model.getLinkProperty(pipeIds, "diameter_m");
demandWeights = max(0.0, model.getNodeProperty(nodeIds, "demand_m3s"));

numberOfPipes = numel(pipeIds);
lowerBounds = 0.01 * ones(1, numberOfPipes);
upperBounds = 1.00 * ones(1, numberOfPipes);
initialDiameters = min(upperBounds, max(lowerBounds, ...
    originalDiameters(:)'));

lastDiameters = [];
lastObjective = [];
lastDetails = [];
[baselineObjective, baselineDetails] = evaluateCached(initialDiameters);
fprintf('Baseline demand-weighted water age: %.3f h\n', ...
    baselineObjective / 3600);
fprintf('Baseline pressure range: %.3f ... %.3f m\n', ...
    min(baselineDetails.pressure_m), max(baselineDetails.pressure_m));
baselineFeasible = isFeasible(baselineDetails);
if ~baselineFeasible
    fprintf(['The initial network violates the 0...60 m pressure constraint; ' ...
             'it is retained as a reference, not as a feasible GA seed.\n']);
end

rng(20260831, 'twister');
if nargin >= 2 && ~isempty(userOptions)
    options = userOptions;
else
    populationSize = max(60, 2 * numberOfPipes);
    options = optimoptions('ga', ...
        'Display', 'iter', ...
        'PopulationSize', populationSize, ...
        'MaxGenerations', 50, ...
        'MaxStallGenerations', 15, ...
        'FunctionTolerance', 1e-5, ...
        'ConstraintTolerance', 1e-4, ...
        'UseParallel', false);
end

% A population drawn uniformly over 0.01...1.00 m is overwhelmingly made
% up of hydraulically infeasible networks. In that situation GA can stop
% after its first generation with exit flag -2. Seed it with several
% feasible, physically related designs while retaining the broad bounds for
% subsequent crossover and mutation.
if isempty(options.InitialPopulationMatrix)
    requestedPopulation = sum(options.PopulationSize);
    targetSeedCount = min(requestedPopulation, max(12, ceil(numberOfPipes / 2)));
    [seedPopulation, seedObjectives] = feasibleInitialPopulation( ...
        targetSeedCount);
    options.InitialPopulationMatrix = seedPopulation;
    fprintf(['Initial GA population: %d feasible seeded designs; ' ...
             'best seeded age %.3f h.\n'], ...
        size(seedPopulation, 1), min(seedObjectives) / 3600);
end

[bestDiameters, ~, exitFlag, output, population, scores] = ga( ...
    @objectiveFunction, numberOfPipes, [], [], [], [], lowerBounds, ...
    upperBounds, @pressureConstraints, options);
[bestObjectiveChecked, bestDetails] = evaluateCached(bestDiameters);
bestObjective = bestObjectiveChecked;

result = struct();
result.network_file = char(networkFile);
result.pipe_ids = pipeIds;
result.node_ids = nodeIds;
result.original_diameters_m = originalDiameters;
result.optimized_diameters_m = bestDiameters(:);
result.baseline_weighted_age_s = baselineObjective;
result.optimized_weighted_age_s = bestObjective;
result.baseline_pressure_m = baselineDetails.pressure_m;
result.optimized_pressure_m = bestDetails.pressure_m;
result.optimized_water_age_s = bestDetails.water_age_s;
result.baseline_feasible = baselineFeasible;
result.optimized_feasible = isFeasible(bestDetails);
result.exit_flag = exitFlag;
result.ga_output = output;
result.final_population = population;
result.final_scores = scores;

outputDirectory = fullfile(sourceDirectory, 'matlab', 'matlab-results');
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
pipeTable = table(pipeIds(:), originalDiameters(:), bestDiameters(:), ...
    'VariableNames', {'pipe_id', 'original_diameter_m', ...
                      'optimized_diameter_m'});
nodeTable = table(nodeIds(:), ...
    baselineDetails.pressure_m(:), bestDetails.pressure_m(:), ...
    bestDetails.water_age_s(:), demandWeights(:), ...
    'VariableNames', {'node_id', 'baseline_pressure_m', ...
                      'optimized_pressure_m', 'optimized_water_age_s', ...
                      'demand_m3s'});
writetable(pipeTable, fullfile(outputDirectory, 'anytown_diameters.csv'));
writetable(nodeTable, fullfile(outputDirectory, 'anytown_nodes.csv'));
save(fullfile(outputDirectory, 'anytown_optimization.mat'), 'result');

figureHandle = figure('Color', 'white', ...
    'Name', 'STACI Anytown diameter optimization');
tiledlayout(figureHandle, 2, 1, 'TileSpacing', 'compact');
nexttile;
diameterBars = [originalDiameters(:), bestDiameters(:)] * 1000;
bar(1:numberOfPipes, diameterBars, 'grouped');
ylabel('Diameter [mm]');
xlabel('Pipe');
xlim([0.25, numberOfPipes + 0.75]);
xticks(1:numberOfPipes);
xticklabels(pipeIds(:));
xtickangle(90);
grid on;
legend({'Initial', 'Optimized'}, 'Location', 'best');
title(sprintf('Weighted age: %.3f h -> %.3f h', ...
    baselineObjective / 3600, bestObjective / 3600));

nexttile;
plot(baselineDetails.pressure_m, 'o-', 'DisplayName', 'Original');
hold on;
plot(bestDetails.pressure_m, '.-', 'DisplayName', 'Optimized');
yline(0, 'r--', '0 m constraint', 'HandleVisibility', 'off');
yline(60, 'r--', '60 m constraint', 'HandleVisibility', 'off');
ylabel('Pressure head [m]');
xlabel('Node index');
grid on;
legend('Location', 'best');
exportgraphics(figureHandle, ...
    fullfile(outputDirectory, 'anytown_optimization.png'), ...
    'Resolution', 180);

fprintf('Optimized demand-weighted water age: %.3f h\n', ...
    bestObjective / 3600);
fprintf('Optimized pressure range: %.3f ... %.3f m\n', ...
    min(bestDetails.pressure_m), max(bestDetails.pressure_m));
fprintf('GA generations: %d; exit flag: %d (%s)\n', ...
    output.generations, exitFlag, output.message);
fprintf('Results: %s\n', outputDirectory);

    function [population, objectives] = feasibleInitialPopulation(targetCount)
        population = zeros(0, numberOfPipes);
        objectives = zeros(0, 1);
        scaleCandidates = unique([1.0, linspace(0.50, 2.50, 41)]);
        for scale = scaleCandidates
            candidate = min(upperBounds, max(lowerBounds, ...
                scale * initialDiameters));
            [candidateObjective, candidateDetails] = evaluateCached(candidate);
            if isFeasible(candidateDetails)
                population(end + 1, :) = candidate; %#ok<AGROW>
                objectives(end + 1, 1) = candidateObjective; %#ok<AGROW>
            end
        end
        if isempty(population)
            error('STACI:NoFeasibleInitialDesign', [ ...
                'No diameter-scaled design satisfies the 0...60 m pressure ' ...
                'constraints. GA was not started.']);
        end

        [~, bestIndex] = min(objectives);
        centre = population(bestIndex, :);
        attempts = 0;
        while size(population, 1) < targetCount && attempts < 20 * targetCount
            attempts = attempts + 1;
            candidate = centre .* exp(0.08 * randn(1, numberOfPipes));
            candidate = min(upperBounds, max(lowerBounds, candidate));
            [candidateObjective, candidateDetails] = evaluateCached(candidate);
            if isFeasible(candidateDetails)
                population(end + 1, :) = candidate; %#ok<AGROW>
                objectives(end + 1, 1) = candidateObjective; %#ok<AGROW>
            end
        end

        [objectives, order] = sort(objectives);
        population = population(order, :);
        keep = min(targetCount, size(population, 1));
        population = population(1:keep, :);
        objectives = objectives(1:keep);
    end

    function feasible = isFeasible(details)
        pressure = details.pressure_m(:);
        feasible = details.converged && all(isfinite(pressure)) && ...
            all(pressure >= -1e-6) && all(pressure <= 60.0 + 1e-6);
    end

    function [objective, details] = evaluateCached(diameters)
        diameters = diameters(:)';
        if isempty(lastDiameters) || ~isequal(diameters, lastDiameters)
            model.setLinkProperty(pipeIds, "diameter_m", diameters);
            hydraulicStatus = model.solveHydraulics();
            lastDetails = struct();
            lastDetails.converged = hydraulicStatus.converged;
            lastDetails.pressure_m = model.getNodeProperty( ...
                nodeIds, "pressure_head_m");
            lastDetails.water_age_s = nan(numel(nodeIds), 1);
            lastDetails.link_flow_m3s = nan(numel(model.Info.link_ids), 1);
            if hydraulicStatus.converged
                nodes = model.nodeResults();
                lastDetails.pressure_m = nodes.pressure_head_m;
                quality = model.solveSteadyWaterAge();
                links = model.linkResults();
                lastDetails.water_age_s = quality.node_water_age_s;
                lastDetails.link_flow_m3s = links.flow_rate_m3s;
                lastObjective = sum(demandWeights .* quality.node_water_age_s) / ...
                    sum(demandWeights);
            else
                lastObjective = 1e12;
            end
            lastDiameters = diameters;
        end
        objective = lastObjective;
        details = lastDetails;
    end

    function objective = objectiveFunction(diameters)
        [objective, details] = evaluateCached(diameters);
        if ~details.converged || ~isfinite(objective)
            objective = 1e12;
        end
    end

    function [inequality, equality] = pressureConstraints(diameters)
        [~, details] = evaluateCached(diameters);
        pressure = details.pressure_m(:);
        if ~details.converged || any(~isfinite(pressure))
            inequality = 1e6 * ones(2 * numel(nodeIds), 1);
        else
            inequality = [-pressure; pressure - 60.0];
        end
        equality = [];
    end
end
