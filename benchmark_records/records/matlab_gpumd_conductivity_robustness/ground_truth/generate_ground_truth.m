% Generate MATLAB ground truth for the synthetic GPUMD robustness record.
% This script reads all 81 public trajectories; it does not import legacy
% Python reference outputs.

scriptPath = mfilename('fullpath');
groundTruthDir = fileparts(scriptPath);
recordRoot = fileparts(groundTruthDir);
inputDir = fullfile(recordRoot, 'inputs');
outputDir = fullfile(groundTruthDir, 'output');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

tracePath = fullfile(outputDir, 'matlab_trace.txt');
if isfile(tracePath)
    delete(tracePath);
end
diary(tracePath);
diaryCleanup = onCleanup(@() diary('off'));

fprintf('GPUMD conductivity ground-truth generation\n');
fprintf('Started: %s\n', char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX')));
fprintf('MATLAB: %s (%s)\n', version, version('-release'));
fprintf('Input directory: %s\n', inputDir);
fprintf('Output directory: %s\n', outputDir);

try
    burnInPs = 600.0;
    robustnessThreshold = 0.05;
    fprintf('SYNTHETIC_DATA=true\n');
    fprintf('BURN_IN_PS=600\n');
    if isempty(ver('stats'))
        error('GPUMD:MissingToolbox', ...
            'Statistics and Machine Learning Toolbox is required.');
    end
    if isempty(ver('optim'))
        error('GPUMD:MissingToolbox', 'Optimization Toolbox is required.');
    end
    normalReferenceCoverage = normcdf(1.96) - normcdf(-1.96);
    fprintf('Statistics Toolbox normal coverage for +/-1.96: %.9f\n', ...
        normalReferenceCoverage);
    fprintf('NORMCDF_COVERAGE=%.15g\n', normalReferenceCoverage);
    metadataPath = fullfile(inputDir, 'runs.csv');
    if ~isfile(metadataPath)
        error('GPUMD:MissingMetadata', 'Missing input file: %s', metadataPath);
    end

    metadata = readtable(metadataPath, 'FileType', 'text', ...
        'Delimiter', ',', 'TextType', 'string');
    requiredMetadata = ["run_id", "material", "size_label", "length_nm", ...
        "drive_label", "driving_force_A_inv", "seed"];
    assertColumns(metadata, requiredMetadata, 'runs.csv');
    if height(metadata) ~= 81
        error('GPUMD:RunCount', 'Expected 81 metadata rows, found %d.', height(metadata));
    end
    if numel(unique(metadata.run_id)) ~= 81
        error('GPUMD:DuplicateRunId', 'runs.csv must contain 81 unique run_id values.');
    end

    nRuns = height(metadata);
    runId = strings(nRuns, 1);
    material = strings(nRuns, 1);
    sizeLabel = strings(nRuns, 1);
    lengthNm = zeros(nRuns, 1);
    driveLabel = strings(nRuns, 1);
    drivingForce = zeros(nRuns, 1);
    seed = zeros(nRuns, 1);
    burnIn = repmat(burnInPs, nRuns, 1);
    nSamples = zeros(nRuns, 1);
    meanTemperature = zeros(nRuns, 1);
    kappaIn = zeros(nRuns, 1);
    kappaOut = zeros(nRuns, 1);
    kappaTotal = zeros(nRuns, 1);
    closureError = zeros(nRuns, 1);

    requiredTrajectory = ["time_ps", "temperature_K", "kappa_in_W_mK", ...
        "kappa_out_W_mK", "kappa_total_W_mK"];

    for i = 1:nRuns
        currentId = metadata.run_id(i);
        trajectoryPath = fullfile(inputDir, 'raw', char(currentId), 'kappa.out');
        if ~isfile(trajectoryPath)
            error('GPUMD:MissingTrajectory', 'Missing trajectory: %s', trajectoryPath);
        end

        samples = readtable(trajectoryPath, 'FileType', 'text', ...
            'Delimiter', ',', 'TextType', 'string');
        assertColumns(samples, requiredTrajectory, char(currentId));

        timePs = samples.time_ps;
        if isempty(timePs) || any(~isfinite(timePs)) || any(diff(timePs) <= 0)
            error('GPUMD:TimeAxis', 'Non-finite or non-monotonic time axis: %s', currentId);
        end

        selected = timePs >= burnInPs;
        if nnz(selected) < 100
            error('GPUMD:PostBurnIn', ...
                'Fewer than 100 samples remain at or after 600 ps: %s', currentId);
        end

        numericBlock = [samples.temperature_K(selected), ...
            samples.kappa_in_W_mK(selected), samples.kappa_out_W_mK(selected), ...
            samples.kappa_total_W_mK(selected)];
        if any(~isfinite(numericBlock), 'all')
            error('GPUMD:NonFinite', 'Non-finite retained values: %s', currentId);
        end

        runId(i) = currentId;
        material(i) = metadata.material(i);
        sizeLabel(i) = metadata.size_label(i);
        lengthNm(i) = metadata.length_nm(i);
        driveLabel(i) = metadata.drive_label(i);
        drivingForce(i) = metadata.driving_force_A_inv(i);
        seed(i) = metadata.seed(i);
        nSamples(i) = nnz(selected);
        meanTemperature(i) = mean(samples.temperature_K(selected));
        kappaIn(i) = mean(samples.kappa_in_W_mK(selected));
        kappaOut(i) = mean(samples.kappa_out_W_mK(selected));
        kappaTotal(i) = mean(samples.kappa_total_W_mK(selected));
        closureError(i) = abs(kappaIn(i) + kappaOut(i) - kappaTotal(i));

        fprintf('Processed %02d/81: %s (%d retained samples)\n', ...
            i, currentId, nSamples(i));
    end

    runEstimates = table(runId, material, sizeLabel, lengthNm, driveLabel, ...
        drivingForce, seed, burnIn, nSamples, meanTemperature, kappaIn, ...
        kappaOut, kappaTotal, closureError, 'VariableNames', { ...
        'run_id', 'material', 'size_label', 'length_nm', 'drive_label', ...
        'driving_force_A_inv', 'seed', 'burn_in_ps', 'n_samples', ...
        'mean_temperature_K', 'kappa_in_W_mK', 'kappa_out_W_mK', ...
        'kappa_total_W_mK', 'component_closure_error_W_mK'});

    runOutputPath = fullfile(outputDir, 'run_estimates.csv');
    writetable(runEstimates, runOutputPath);
    fprintf('TRAJECTORIES_PROCESSED=81\n');
    fprintf('POST_BURN_IN_UNIT=independent_run\n');

    materialNames = ["CsSnBr3", "Cs2SnBr6", "hetero_N2"];
    materialsPayload = struct();
    sizePayload = struct();
    drivePayload = struct();
    sizeConvergencePayload = struct();
    materialTotalMeans = zeros(numel(materialNames), 1);
    fitModel = @(parameters, lengthValues) ...
        parameters(1) - parameters(2) ./ lengthValues;
    fitOptions = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'Algorithm', 'trust-region-reflective', ...
        'FunctionTolerance', 1.0e-14, ...
        'StepTolerance', 1.0e-14, ...
        'OptimalityTolerance', 1.0e-14, ...
        'MaxIterations', 1000, ...
        'MaxFunctionEvaluations', 5000);

    for i = 1:numel(materialNames)
        name = materialNames(i);
        mask = runEstimates.material == name;
        if nnz(mask) ~= 27
            error('GPUMD:MaterialCount', ...
                'Expected 27 runs for %s, found %d.', name, nnz(mask));
        end

        totals = runEstimates.kappa_total_W_mK(mask);
        totalMean = mean(totals);
        standardError = std(totals, 0) / sqrt(numel(totals));
        materialTotalMeans(i) = totalMean;
        key = char(name);

        materialsPayload.(key) = struct( ...
            'kappa_in_W_mK', mean(runEstimates.kappa_in_W_mK(mask)), ...
            'kappa_out_W_mK', mean(runEstimates.kappa_out_W_mK(mask)), ...
            'kappa_total_W_mK', totalMean, ...
            'standard_error_W_mK', standardError, ...
            'ci95_low_W_mK', totalMean - 1.96 * standardError, ...
            'ci95_high_W_mK', totalMean + 1.96 * standardError, ...
            'n_runs', nnz(mask));

        sizePayload.(key) = computeRobustness(runEstimates(mask, :), ...
            'size_label', ["small", "medium", "large"], robustnessThreshold);
        drivePayload.(key) = computeRobustness(runEstimates(mask, :), ...
            'drive_label', ["low", "medium", "high"], robustnessThreshold);

        lengthsNm = [8.0; 16.0; 23.0];
        observedMeans = [ ...
            sizePayload.(key).group_means_W_mK.small; ...
            sizePayload.(key).group_means_W_mK.medium; ...
            sizePayload.(key).group_means_W_mK.large];
        spread = max(observedMeans) - min(observedMeans);
        kInfInitial = max(observedMeans) + max(spread, 0.01 * mean(observedMeans));
        aInitial = max(mean((kInfInitial - observedMeans) .* lengthsNm), eps);
        initialParameters = [kInfInitial, aInitial];
        lowerBounds = [0.0, 0.0];
        upperBounds = [Inf, Inf];

        [fitParameters, residualNorm, ~, exitFlag] = lsqcurvefit( ...
            fitModel, initialParameters, lengthsNm, observedMeans, ...
            lowerBounds, upperBounds, fitOptions);
        fittedMeans = fitModel(fitParameters, lengthsNm);
        fitRmse = sqrt(mean((fittedMeans - observedMeans) .^ 2));
        relativeRmse = fitRmse / mean(observedMeans);
        fitConverged = exitFlag > 0 && ...
            fitParameters(1) >= max(observedMeans) - 1.0e-12 && ...
            fitParameters(2) >= 0 && relativeRmse <= 0.01;

        sizeConvergencePayload.(key) = struct( ...
            'model', 'kappa(L)=k_inf-a/L', ...
            'lengths_nm', reshape(lengthsNm, 1, []), ...
            'observed_group_means_W_mK', reshape(observedMeans, 1, []), ...
            'fitted_group_means_W_mK', reshape(fittedMeans, 1, []), ...
            'k_inf_W_mK', fitParameters(1), ...
            'a_W_nm_mK', fitParameters(2), ...
            'rmse_W_mK', fitRmse, ...
            'relative_rmse', relativeRmse, ...
            'residual_norm', residualNorm, ...
            'exitflag', exitFlag, ...
            'converged', fitConverged);
        fprintf(['Size fit %s: k_inf=%.9g, a=%.9g, RMSE=%.9g, ', ...
            'relative_RMSE=%.9g, exitflag=%d\n'], key, fitParameters(1), ...
            fitParameters(2), fitRmse, relativeRmse, exitFlag);
    end
    fprintf('LSQCURVEFIT_MATERIALS=3\n');

    [~, orderIndex] = sort(materialTotalMeans, 'descend');
    physicalOrder = cellstr(materialNames(orderIndex));
    expectedOrder = ["CsSnBr3"; "Cs2SnBr6"; "hetero_N2"];
    sizeFlags = arrayfun(@(i) sizePayload.(char(materialNames(i))).robust_within_5_percent, ...
        1:numel(materialNames));
    driveFlags = arrayfun(@(i) drivePayload.(char(materialNames(i))).robust_within_5_percent, ...
        1:numel(materialNames));

    summary = struct();
    summary.schema_version = '1.0';
    summary.synthetic_data = true;
    summary.method = struct( ...
        'burn_in_ps', burnInPs, ...
        'run_estimator', 'arithmetic mean for time_ps >= 600 ps', ...
        'uncertainty', 'sample standard deviation across 27 runs divided by sqrt(27)', ...
        'confidence_interval', 'mean +/- 1.96 * standard_error', ...
        'normal_reference_coverage', normalReferenceCoverage, ...
        'robustness_threshold', robustnessThreshold, ...
        'ensemble', 'all sizes, drives, and independent seeds');
    summary.materials = materialsPayload;
    summary.physical_order = physicalOrder;
    summary.size_robustness = sizePayload;
    summary.drive_robustness = drivePayload;
    summary.size_convergence_fit = sizeConvergencePayload;
    summary.conclusions = struct( ...
        'size_invariant_within_5_percent', all(sizeFlags), ...
        'drive_invariant_within_5_percent', all(driveFlags), ...
        'ordering_matches_expected_group_order', ...
        isequal(string(physicalOrder(:)), expectedOrder));

    jsonText = jsonencode(summary, 'PrettyPrint', true);
    summaryPath = fullfile(outputDir, 'conductivity_summary.json');
    fileId = fopen(summaryPath, 'w', 'n', 'UTF-8');
    if fileId < 0
        error('GPUMD:WriteJSON', 'Cannot open output file: %s', summaryPath);
    end
    fileCleanup = onCleanup(@() fclose(fileId));
    fprintf(fileId, '%s\n', jsonText);
    clear fileCleanup;

    save(fullfile(outputDir, 'ground_truth.mat'), 'summary', 'runEstimates', ...
        'materialNames', 'burnInPs', 'robustnessThreshold', ...
        'sizeConvergencePayload', '-v7');

    fprintf('Wrote: %s\n', summaryPath);
    fprintf('Wrote: %s\n', runOutputPath);
    fprintf('Wrote: %s\n', fullfile(outputDir, 'ground_truth.mat'));
    fprintf('OUTPUT_JSON=conductivity_summary.json\n');
    fprintf('OUTPUT_CSV=run_estimates.csv\n');
    fprintf('Completed successfully: %s\n', char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX')));
catch ME
    fprintf(2, 'GROUND-TRUTH GENERATION FAILED: %s\n', ME.message);
    for i = 1:numel(ME.stack)
        fprintf(2, '  at %s line %d\n', ME.stack(i).name, ME.stack(i).line);
    end
    clear diaryCleanup;
    rethrow(ME);
end

clear diaryCleanup;

function assertColumns(inputTable, requiredNames, sourceName)
actualNames = string(inputTable.Properties.VariableNames);
missing = requiredNames(~ismember(requiredNames, actualNames));
if ~isempty(missing)
    error('GPUMD:MissingColumn', 'Missing columns in %s: %s', ...
        sourceName, strjoin(missing, ', '));
end
end

function result = computeRobustness(runTable, groupVariable, labels, threshold)
groupMeans = struct();
values = zeros(numel(labels), 1);
groups = runTable.(groupVariable);

for i = 1:numel(labels)
    mask = groups == labels(i);
    if ~any(mask)
        error('GPUMD:MissingGroup', 'Missing group %s for %s.', labels(i), groupVariable);
    end
    values(i) = mean(runTable.kappa_total_W_mK(mask));
    groupMeans.(char(labels(i))) = values(i);
end

center = mean(values);
maximumDeviation = max(abs(values - center) ./ abs(center));
result = struct( ...
    'group_means_W_mK', groupMeans, ...
    'max_relative_deviation', maximumDeviation, ...
    'robust_within_5_percent', maximumDeviation <= threshold);
end
