function solution(input_dir, submission_dir)
%SOLUTION Analyze independent GPUMD-style conductivity trajectories.

burnInPs = 600;
requiredRunColumns = [ ...
    "run_id", "material", "size_label", "atom_count", "length_nm", ...
    "drive_label", "driving_force_A_inv", "seed", ...
    "temperature_target_K", "sample_interval_ps", "total_time_ps"];
requiredTrajectoryColumns = [ ...
    "time_ps", "temperature_K", "kappa_in_W_mK", ...
    "kappa_out_W_mK", "kappa_total_W_mK"];

runsPath = fullfile(input_dir, 'runs.csv');
if ~isfile(runsPath)
    error('solution:MissingRunsFile', 'Missing input file: %s', runsPath);
end

runs = readtable(runsPath, 'Delimiter', ',');
require_columns(runs, requiredRunColumns, 'runs.csv');
if height(runs) ~= 81
    error('solution:UnexpectedRunCount', ...
        'Expected 81 rows in runs.csv, found %d.', height(runs));
end

runIds = as_text(runs.run_id, 'run_id');
materials = as_text(runs.material, 'material');
sizeLabels = as_text(runs.size_label, 'size_label');
lengthsNm = as_numeric(runs.length_nm, 'length_nm');
driveLabels = as_text(runs.drive_label, 'drive_label');
drivingForces = as_numeric(runs.driving_force_A_inv, 'driving_force_A_inv');
seeds = as_numeric(runs.seed, 'seed');

% Validate the remaining declared run metadata even though it is not emitted.
as_numeric(runs.atom_count, 'atom_count');
as_numeric(runs.temperature_target_K, 'temperature_target_K');
as_numeric(runs.sample_interval_ps, 'sample_interval_ps');
as_numeric(runs.total_time_ps, 'total_time_ps');

nRuns = height(runs);
if numel(unique(runIds)) ~= nRuns
    error('solution:DuplicateRunId', 'runs.csv contains duplicate run_id values.');
end
if any(seeds ~= round(seeds))
    error('solution:InvalidSeed', 'All seed values must be integers.');
end

nSamples = zeros(nRuns, 1);
meanTemperature = zeros(nRuns, 1);
meanKappaIn = zeros(nRuns, 1);
meanKappaOut = zeros(nRuns, 1);
meanKappaTotal = zeros(nRuns, 1);

for i = 1:nRuns
    trajectoryPath = fullfile(input_dir, 'raw', char(runIds(i)), 'kappa.out');
    if ~isfile(trajectoryPath)
        error('solution:MissingTrajectory', ...
            'Missing trajectory for run_id %s.', runIds(i));
    end

    trajectory = readtable(trajectoryPath, 'Delimiter', ',');
    require_columns(trajectory, requiredTrajectoryColumns, trajectoryPath);

    timePs = as_numeric(trajectory.time_ps, 'time_ps');
    temperatureK = as_numeric(trajectory.temperature_K, 'temperature_K');
    kappaIn = as_numeric(trajectory.kappa_in_W_mK, 'kappa_in_W_mK');
    kappaOut = as_numeric(trajectory.kappa_out_W_mK, 'kappa_out_W_mK');
    kappaTotal = as_numeric(trajectory.kappa_total_W_mK, 'kappa_total_W_mK');

    rowCount = numel(timePs);
    if any([numel(temperatureK), numel(kappaIn), numel(kappaOut), ...
            numel(kappaTotal)] ~= rowCount)
        error('solution:TrajectoryLengthMismatch', ...
            'Trajectory columns have inconsistent lengths for run_id %s.', runIds(i));
    end
    if rowCount < 2 || any(diff(timePs) <= 0)
        error('solution:NonIncreasingTime', ...
            'time_ps must be strictly increasing for run_id %s.', runIds(i));
    end

    retained = timePs >= burnInPs;
    nSamples(i) = nnz(retained);
    if nSamples(i) < 100
        error('solution:TooFewPostBurnInSamples', ...
            'run_id %s has only %d samples at or after 600 ps.', ...
            runIds(i), nSamples(i));
    end

    meanTemperature(i) = mean(temperatureK(retained));
    meanKappaIn(i) = mean(kappaIn(retained));
    meanKappaOut(i) = mean(kappaOut(retained));
    meanKappaTotal(i) = mean(kappaTotal(retained));
end

closureError = abs(meanKappaIn + meanKappaOut - meanKappaTotal);
runEstimates = table( ...
    runIds, materials, sizeLabels, lengthsNm, driveLabels, drivingForces, ...
    seeds, repmat(burnInPs, nRuns, 1), nSamples, meanTemperature, ...
    meanKappaIn, meanKappaOut, meanKappaTotal, closureError, ...
    'VariableNames', { ...
        'run_id', 'material', 'size_label', 'length_nm', 'drive_label', ...
        'driving_force_A_inv', 'seed', 'burn_in_ps', 'n_samples', ...
        'mean_temperature_K', 'kappa_in_W_mK', 'kappa_out_W_mK', ...
        'kappa_total_W_mK', 'component_closure_error_W_mK'});

materialNames = ["CsSnBr3", "Cs2SnBr6", "hetero_N2"];
sizeLevels = ["small", "medium", "large"];
driveLevels = ["low", "medium", "high"];
fitLengthsNm = [8; 16; 23];

validate_factor(materials, materialNames, 27, 'material');
validate_factor(sizeLabels, sizeLevels, 27, 'size_label');
validate_factor(driveLabels, driveLevels, 27, 'drive_label');
for j = 1:numel(sizeLevels)
    sizeMask = sizeLabels == sizeLevels(j);
    if any(abs(lengthsNm(sizeMask) - fitLengthsNm(j)) > 1e-12)
        error('solution:SizeLengthMismatch', ...
            'size_label %s does not map uniquely to length %.15g nm.', ...
            sizeLevels(j), fitLengthsNm(j));
    end
end

% This call is intentionally evaluated at runtime using Statistics Toolbox.
normalCoverage = normcdf(1.96) - normcdf(-1.96);
if ~isfinite(normalCoverage)
    error('solution:InvalidNormalCoverage', 'normcdf returned a nonfinite value.');
end

summary = struct();
summary.schema_version = '1.0';
summary.synthetic_data = true;
summary.method = struct( ...
    'burn_in_ps', burnInPs, ...
    'run_estimator', 'arithmetic_mean_post_burn_in', ...
    'ensemble', 'independent_run', ...
    'normal_reference_coverage', normalCoverage);
summary.materials = struct();
summary.size_robustness = struct();
summary.drive_robustness = struct();
summary.size_convergence_fit = struct();

materialTotalMeans = zeros(1, numel(materialNames));
sizeRobustFlags = false(1, numel(materialNames));
driveRobustFlags = false(1, numel(materialNames));

for i = 1:numel(materialNames)
    materialName = materialNames(i);
    materialMask = materials == materialName;
    if nnz(materialMask) ~= 27
        error('solution:MaterialRunCount', ...
            'Material %s must contain exactly 27 independent runs.', materialName);
    end

    inValues = meanKappaIn(materialMask);
    outValues = meanKappaOut(materialMask);
    totalValues = meanKappaTotal(materialMask);
    totalMean = mean(totalValues);
    standardError = std(totalValues, 0) / sqrt(numel(totalValues));

    materialSummary = struct();
    materialSummary.kappa_in_W_mK = mean(inValues);
    materialSummary.kappa_out_W_mK = mean(outValues);
    materialSummary.kappa_total_W_mK = totalMean;
    materialSummary.standard_error_W_mK = standardError;
    materialSummary.ci95_low_W_mK = totalMean - 1.96 * standardError;
    materialSummary.ci95_high_W_mK = totalMean + 1.96 * standardError;
    materialSummary.n_runs = numel(totalValues);
    summary.materials.(char(materialName)) = materialSummary;
    materialTotalMeans(i) = totalMean;

    [sizeSummary, sizeMeans, sizeRobust] = robustness_summary( ...
        meanKappaTotal, materialMask, sizeLabels, sizeLevels, 9);
    [driveSummary, ~, driveRobust] = robustness_summary( ...
        meanKappaTotal, materialMask, driveLabels, driveLevels, 9);
    summary.size_robustness.(char(materialName)) = sizeSummary;
    summary.drive_robustness.(char(materialName)) = driveSummary;
    sizeRobustFlags(i) = sizeRobust;
    driveRobustFlags(i) = driveRobust;

    fitSummary = fit_size_convergence(fitLengthsNm, sizeMeans);
    summary.size_convergence_fit.(char(materialName)) = fitSummary;
end

[~, orderIndices] = sort(materialTotalMeans, 'descend');
physicalOrder = materialNames(orderIndices);
summary.physical_order = cellstr(physicalOrder);
expectedOrder = ["CsSnBr3", "Cs2SnBr6", "hetero_N2"];

summary.conclusions = struct( ...
    'size_invariant_within_5_percent', logical(all(sizeRobustFlags)), ...
    'drive_invariant_within_5_percent', logical(all(driveRobustFlags)), ...
    'ordering_matches_expected_group_order', logical(isequal(physicalOrder, expectedOrder)));

assert_finite_numeric_content(summary, 'conductivity_summary');
jsonText = jsonencode(summary);

if ~isfolder(submission_dir)
    [created, message] = mkdir(submission_dir);
    if ~created
        error('solution:CannotCreateSubmissionDirectory', ...
            'Cannot create submission directory: %s', message);
    end
end

csvName = 'run_estimates.csv';
jsonName = 'conductivity_summary.json';
traceName = 'matlab_trace.txt';
writetable(runEstimates, fullfile(submission_dir, csvName), 'Delimiter', ',');
write_text_file(fullfile(submission_dir, jsonName), [jsonText, newline]);

traceText = sprintf([ ...
    'MATLAB_ANALYSIS_STATUS=completed\n' ...
    'VALIDATION_STATUS=passed\n' ...
    'SYNTHETIC_DATA=true\n' ...
    'TRAJECTORIES_PROCESSED=%d\n' ...
    'BURN_IN_PS=%.0f\n' ...
    'POST_BURN_IN_UNIT=independent_run\n' ...
    'NORMCDF_COVERAGE=%.17g\n' ...
    'LSQCURVEFIT_MATERIALS=%d\n' ...
    'OUTPUT_JSON=%s\n' ...
    'OUTPUT_CSV=%s\n' ...
    'OUTPUT_STATUS=written\n'], ...
    nRuns, burnInPs, normalCoverage, numel(materialNames), jsonName, csvName);
write_text_file(fullfile(submission_dir, traceName), traceText);
end


function require_columns(inputTable, requiredColumns, sourceName)
availableColumns = string(inputTable.Properties.VariableNames);
missingColumns = requiredColumns(~ismember(requiredColumns, availableColumns));
if ~isempty(missingColumns)
    error('solution:MissingColumns', '%s is missing columns: %s', ...
        sourceName, strjoin(missingColumns, ', '));
end
end


function values = as_text(column, columnName)
values = strip(string(column));
values = values(:);
if any(ismissing(values) | strlength(values) == 0)
    error('solution:InvalidTextColumn', ...
        'Column %s contains a missing or empty value.', columnName);
end
end


function values = as_numeric(column, columnName)
if isnumeric(column) || islogical(column)
    values = double(column);
else
    values = str2double(string(column));
end
values = values(:);
if any(~isfinite(values))
    error('solution:InvalidNumericColumn', ...
        'Column %s contains a nonfinite or nonnumeric value.', columnName);
end
end


function validate_factor(observed, expectedLevels, expectedCount, factorName)
unexpected = setdiff(unique(observed), expectedLevels);
if ~isempty(unexpected)
    error('solution:UnexpectedFactorLevel', ...
        '%s contains unexpected levels: %s', factorName, strjoin(unexpected, ', '));
end
for i = 1:numel(expectedLevels)
    count = nnz(observed == expectedLevels(i));
    if count ~= expectedCount
        error('solution:FactorCount', ...
            '%s level %s has %d rows; expected %d.', ...
            factorName, expectedLevels(i), count, expectedCount);
    end
end
end


function [result, groupMeans, robust] = robustness_summary( ...
        response, materialMask, observedLevels, expectedLevels, expectedCount)
groupMeans = zeros(numel(expectedLevels), 1);
groupMeanFields = struct();
for i = 1:numel(expectedLevels)
    groupMask = materialMask & observedLevels == expectedLevels(i);
    if nnz(groupMask) ~= expectedCount
        error('solution:WithinMaterialFactorCount', ...
            'Level %s has %d runs within a material; expected %d.', ...
            expectedLevels(i), nnz(groupMask), expectedCount);
    end
    groupMeans(i) = mean(response(groupMask));
    groupMeanFields.(char(expectedLevels(i))) = groupMeans(i);
end

center = mean(groupMeans);
if ~isfinite(center) || center == 0
    error('solution:InvalidRobustnessCenter', ...
        'The arithmetic center for robustness must be finite and nonzero.');
end
maxRelativeDeviation = max(abs(groupMeans - center) / center);
robust = maxRelativeDeviation <= 0.05;

result = struct();
result.group_means_W_mK = groupMeanFields;
result.max_relative_deviation = maxRelativeDeviation;
result.robust_within_5_percent = logical(robust);
end


function result = fit_size_convergence(lengthsNm, observedMeans)
lengthsNm = lengthsNm(:);
observedMeans = observedMeans(:);
if numel(lengthsNm) ~= 3 || numel(observedMeans) ~= 3
    error('solution:FitPointCount', ...
        'Size convergence fitting requires exactly three size points.');
end
if any(lengthsNm <= 0) || any(~isfinite(observedMeans)) || mean(observedMeans) == 0
    error('solution:InvalidFitData', 'Size convergence fit data are invalid.');
end

kInf0 = max(observedMeans) + max( ...
    max(observedMeans) - min(observedMeans), 0.01 * mean(observedMeans));
a0 = max(mean((kInf0 - observedMeans) .* lengthsNm), eps);
modelFunction = @(p, L) p(1) - p(2) ./ L;
options = optimoptions('lsqcurvefit', 'Display', 'off');
[parameters, ~, ~, exitflag] = lsqcurvefit( ...
    modelFunction, [kInf0, a0], lengthsNm, observedMeans, ...
    [0, 0], [Inf, Inf], options);

fittedMeans = modelFunction(parameters, lengthsNm);
rmse = sqrt(mean((fittedMeans - observedMeans) .^ 2));
relativeRmse = rmse / mean(observedMeans);
if any(~isfinite([parameters(:); fittedMeans(:); rmse; relativeRmse; exitflag]))
    error('solution:NonfiniteFit', 'lsqcurvefit returned a nonfinite result.');
end
converged = exitflag > 0 && ...
    parameters(1) >= max(observedMeans) && ...
    parameters(2) >= 0 && relativeRmse <= 0.01;

result = struct();
result.model = 'kappa(L)=k_inf-a/L';
result.lengths_nm = lengthsNm.';
result.observed_group_means_W_mK = observedMeans.';
result.fitted_group_means_W_mK = fittedMeans.';
result.k_inf_W_mK = parameters(1);
result.a_W_nm_mK = parameters(2);
result.rmse_W_mK = rmse;
result.relative_rmse = relativeRmse;
result.exitflag = double(exitflag);
result.converged = logical(converged);
end


function assert_finite_numeric_content(value, valueName)
if isnumeric(value)
    if any(~isfinite(value(:)))
        error('solution:NonfiniteOutput', ...
            '%s contains a nonfinite numeric value.', valueName);
    end
elseif isstruct(value)
    fields = fieldnames(value);
    for elementIndex = 1:numel(value)
        for fieldIndex = 1:numel(fields)
            childName = sprintf('%s.%s', valueName, fields{fieldIndex});
            assert_finite_numeric_content( ...
                value(elementIndex).(fields{fieldIndex}), childName);
        end
    end
elseif iscell(value)
    for i = 1:numel(value)
        assert_finite_numeric_content(value{i}, sprintf('%s{%d}', valueName, i));
    end
end
end


function write_text_file(filePath, textValue)
fileId = fopen(filePath, 'w');
if fileId < 0
    error('solution:CannotWriteFile', 'Cannot open output file: %s', filePath);
end
cleanup = onCleanup(@() fclose(fileId)); %#ok<NASGU>
count = fwrite(fileId, char(textValue), 'char');
if count ~= numel(char(textValue))
    error('solution:ShortWrite', 'Incomplete write to output file: %s', filePath);
end
end
