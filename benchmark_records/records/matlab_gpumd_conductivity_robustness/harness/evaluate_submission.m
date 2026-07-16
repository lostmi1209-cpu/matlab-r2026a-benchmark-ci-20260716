function report = evaluate_submission(submission_dir, ground_truth_dir, report_path)
%EVALUATE_SUBMISSION Deterministically grade the GPUMD robustness outputs.
%   REPORT = EVALUATE_SUBMISSION(SUBMISSION_DIR, GROUND_TRUTH_DIR, REPORT_PATH)
%   returns a report struct and writes the same report as JSON.

arguments
    submission_dir (1, 1) string
    ground_truth_dir (1, 1) string
    report_path (1, 1) string
end

maxScore = 100.0;
score = 0.0;
details = struct('id', {}, 'points', {}, 'earned', {}, 'message', {});

candidateSummaryPath = fullfile(submission_dir, 'conductivity_summary.json');
candidateRunsPath = fullfile(submission_dir, 'run_estimates.csv');
candidateTracePath = fullfile(submission_dir, 'matlab_trace.txt');
truthSummaryPath = fullfile(ground_truth_dir, 'conductivity_summary.json');
truthRunsPath = fullfile(ground_truth_dir, 'run_estimates.csv');

if ~isfile(truthSummaryPath) || ~isfile(truthRunsPath)
    error('GPUMD:MissingGroundTruth', ...
        'Ground Truth files are missing from %s.', ground_truth_dir);
end

truthSummary = jsondecode(fileread(truthSummaryPath));
truthRuns = readtable(truthRunsPath, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');

summaryExists = isfile(candidateSummaryPath);
runsExist = isfile(candidateRunsPath);
traceExists = isfile(candidateTracePath);
summaryParsed = false;
runsParsed = false;
candidateSummary = struct();
candidateRuns = table();
traceValid = false;

if summaryExists
    try
        candidateSummary = jsondecode(fileread(candidateSummaryPath));
        summaryParsed = true;
    catch
        summaryParsed = false;
    end
end

if runsExist
    try
        candidateRuns = readtable(candidateRunsPath, 'FileType', 'text', ...
            'Delimiter', ',', 'TextType', 'string');
        runsParsed = true;
    catch
        runsParsed = false;
    end
end

if traceExists
    try
        traceText = string(fileread(candidateTracePath));
        traceLines = strtrim(splitlines(traceText));
        requiredTraceLines = [ ...
            "SYNTHETIC_DATA=true", ...
            "TRAJECTORIES_PROCESSED=81", ...
            "BURN_IN_PS=600", ...
            "POST_BURN_IN_UNIT=independent_run", ...
            "LSQCURVEFIT_MATERIALS=3", ...
            "OUTPUT_JSON=conductivity_summary.json", ...
            "OUTPUT_CSV=run_estimates.csv"];
        normLineMask = startsWith(traceLines, "NORMCDF_COVERAGE=");
        normLineValid = false;
        if nnz(normLineMask) == 1
            traceCoverage = str2double(extractAfter(traceLines(normLineMask), ...
                "NORMCDF_COVERAGE="));
            normLineValid = scalarMatch(traceCoverage, ...
                truthSummary.method.normal_reference_coverage, 1.0e-12);
        end
        traceValid = strlength(strtrim(traceText)) > 0 && ...
            all(ismember(requiredTraceLines, traceLines)) && normLineValid;
    catch
        traceValid = false;
    end
end

earned = double(summaryExists) + double(runsExist) + double(traceExists) + ...
    double(summaryParsed) + double(runsParsed);
schemaOk = false;
if summaryParsed && isstruct(candidateSummary) && ...
        isfield(candidateSummary, 'schema_version')
    schemaValue = string(candidateSummary.schema_version);
    schemaOk = isscalar(schemaValue) && schemaValue == "1.0";
end
earned = earned + double(schemaOk);
syntheticOk = false;
if summaryParsed && isstruct(candidateSummary) && ...
        isfield(candidateSummary, 'synthetic_data') && ...
        islogical(candidateSummary.synthetic_data) && ...
        isscalar(candidateSummary.synthetic_data)
    syntheticOk = candidateSummary.synthetic_data;
end
earned = earned + double(syntheticOk);
earned = earned + double(traceValid);
[score, details] = award(score, details, 'files_and_schema', 8.0, earned, ...
    'Required JSON/CSV/Trace, parsing, schema fields, and trace markers.');

matchedTruthRows = false(height(truthRuns), 1);
candidateLocations = zeros(height(truthRuns), 1);
candidateIdsUnique = false;
requiredRunColumns = ["run_id", "burn_in_ps", "n_samples", ...
    "mean_temperature_K", "kappa_in_W_mK", "kappa_out_W_mK", ...
    "kappa_total_W_mK", "component_closure_error_W_mK"];

if runsParsed && hasColumns(candidateRuns, requiredRunColumns)
    candidateRunIds = string(candidateRuns.run_id);
    truthRunIds = string(truthRuns.run_id);
    candidateIdsUnique = numel(unique(candidateRunIds)) == height(candidateRuns);
    [matchedTruthRows, candidateLocations] = ismember(truthRunIds, candidateRunIds);
    overlap = nnz(matchedTruthRows) / height(truthRuns);
    exact = height(candidateRuns) == height(truthRuns) && candidateIdsUnique && all(matchedTruthRows);
    earned = 6.0 * overlap + 2.0 * double(exact);
else
    earned = 0.0;
end
[score, details] = award(score, details, 'run_coverage', 8.0, earned, ...
    'Coverage of 81 unique run_id values.');

if runsParsed && candidateIdsUnique && hasColumns(candidateRuns, requiredRunColumns)
    validTruth = find(matchedTruthRows);
    validCandidate = candidateLocations(matchedTruthRows);
    burnMatches = numericMatches(candidateRuns.burn_in_ps(validCandidate), ...
        truthRuns.burn_in_ps(validTruth), 1.0e-9);
    sampleMatches = numericMatches(candidateRuns.n_samples(validCandidate), ...
        truthRuns.n_samples(validTruth), 0.0);
    denominator = height(truthRuns);
    earned = 4.0 * nnz(burnMatches) / denominator + ...
        4.0 * nnz(sampleMatches) / denominator;
else
    earned = 0.0;
end
[score, details] = award(score, details, 'equilibration_rule', 8.0, earned, ...
    '600 ps burn-in and retained sample counts.');

if runsParsed && candidateIdsUnique && hasColumns(candidateRuns, requiredRunColumns)
    validTruth = find(matchedTruthRows);
    validCandidate = candidateLocations(matchedTruthRows);
    numericFields = ["mean_temperature_K", "kappa_in_W_mK", ...
        "kappa_out_W_mK", "kappa_total_W_mK", ...
        "component_closure_error_W_mK"];
    tolerances = [1.0e-3, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6];
    earned = 0.0;
    for i = 1:numel(numericFields)
        field = numericFields(i);
        candidateValues = candidateRuns.(field);
        truthValues = truthRuns.(field);
        matches = numericMatches(candidateValues(validCandidate), ...
            truthValues(validTruth), tolerances(i));
        earned = earned + 3.6 * nnz(matches) / height(truthRuns);
    end
else
    earned = 0.0;
end
[score, details] = award(score, details, 'run_level_estimates', 18.0, earned, ...
    'Run-level temperature, conductivity components, total, and closure.');

materialNames = ["CsSnBr3", "Cs2SnBr6", "hetero_N2"];
materialFields = ["kappa_in_W_mK", "kappa_out_W_mK", ...
    "kappa_total_W_mK", "standard_error_W_mK", ...
    "ci95_low_W_mK", "ci95_high_W_mK", "n_runs"];
correctMaterialFields = 0;
totalMaterialFields = numel(materialNames) * numel(materialFields);

if summaryParsed && isstruct(candidateSummary) && ...
        isfield(candidateSummary, 'materials') && isstruct(candidateSummary.materials)
    for i = 1:numel(materialNames)
        material = char(materialNames(i));
        if ~isfield(candidateSummary.materials, material) || ...
                ~isfield(truthSummary.materials, material)
            continue;
        end
        candidateMaterial = candidateSummary.materials.(material);
        truthMaterial = truthSummary.materials.(material);
        if ~isstruct(candidateMaterial)
            continue;
        end
        for j = 1:numel(materialFields)
            field = char(materialFields(j));
            if ~isfield(candidateMaterial, field) || ~isfield(truthMaterial, field)
                continue;
            end
            tolerance = 1.0e-6;
            if strcmp(field, 'n_runs')
                tolerance = 0.0;
            end
            if scalarMatch(candidateMaterial.(field), truthMaterial.(field), tolerance)
                correctMaterialFields = correctMaterialFields + 1;
            end
        end
    end
end
earned = 18.0 * correctMaterialFields / totalMaterialFields;
[score, details] = award(score, details, 'material_aggregation', 18.0, earned, ...
    'Material means, uncertainty, confidence intervals, and run counts.');

closureEarned = 0.0;
if summaryParsed && isstruct(candidateSummary) && ...
        isfield(candidateSummary, 'materials') && isstruct(candidateSummary.materials)
    for i = 1:numel(materialNames)
        material = char(materialNames(i));
        if ~isfield(candidateSummary.materials, material)
            continue;
        end
        entry = candidateSummary.materials.(material);
        if ~isstruct(entry)
            continue;
        end
        required = {'kappa_in_W_mK', 'kappa_out_W_mK', 'kappa_total_W_mK'};
        if all(isfield(entry, required))
            values = [double(entry.kappa_in_W_mK), double(entry.kappa_out_W_mK), ...
                double(entry.kappa_total_W_mK)];
            if all(isfinite(values)) && abs(values(1) + values(2) - values(3)) <= 0.005
                closureEarned = closureEarned + 2.0 / 3.0;
            end
        end
    end
end

orderEarned = 0.0;
if summaryParsed && isstruct(candidateSummary) && ...
        isfield(candidateSummary, 'physical_order') && ...
        isfield(truthSummary, 'physical_order')
    try
        if isequal(string(candidateSummary.physical_order(:)), ...
                string(truthSummary.physical_order(:)))
            orderEarned = 2.0;
        end
    catch
        orderEarned = 0.0;
    end
end

conclusionEarned = 0.0;
conclusionFields = ["size_invariant_within_5_percent", ...
    "drive_invariant_within_5_percent", "ordering_matches_expected_group_order"];
if summaryParsed && isstruct(candidateSummary) && ...
        isfield(candidateSummary, 'conclusions') && ...
        isstruct(candidateSummary.conclusions) && ...
        isfield(truthSummary, 'conclusions')
    for i = 1:numel(conclusionFields)
        field = char(conclusionFields(i));
        if isfield(candidateSummary.conclusions, field) && ...
                isfield(truthSummary.conclusions, field) && ...
                isequal(candidateSummary.conclusions.(field), truthSummary.conclusions.(field))
            conclusionEarned = conclusionEarned + 1.0;
        end
    end
end
coverageEarned = 0.0;
if summaryParsed && isstruct(candidateSummary) && ...
        isfield(candidateSummary, 'method') && isstruct(candidateSummary.method) && ...
        isfield(candidateSummary.method, 'normal_reference_coverage') && ...
        isfield(truthSummary, 'method') && ...
        isfield(truthSummary.method, 'normal_reference_coverage') && ...
        scalarMatch(candidateSummary.method.normal_reference_coverage, ...
        truthSummary.method.normal_reference_coverage, 1.0e-12)
    coverageEarned = 1.0;
end
earned = closureEarned + orderEarned + conclusionEarned + coverageEarned;
[score, details] = award(score, details, 'statistical_consistency', 8.0, earned, ...
    'Arithmetic closure, mean ordering, normcdf coverage, and conclusion flags.');

earned = robustnessScore(candidateSummary, truthSummary, ...
    'size_robustness', materialNames, ["small", "medium", "large"]);
[score, details] = award(score, details, 'size_robustness', 10.0, earned, ...
    'Size-group means, maximum relative deviations, and 5 percent flags.');

earned = robustnessScore(candidateSummary, truthSummary, ...
    'drive_robustness', materialNames, ["low", "medium", "high"]);
[score, details] = award(score, details, 'drive_robustness', 10.0, earned, ...
    'Drive-group means, maximum relative deviations, and 5 percent flags.');

earned = sizeConvergenceScore(candidateSummary, truthSummary, materialNames);
[score, details] = award(score, details, 'size_convergence_optimization', ...
    12.0, earned, ...
    'Constrained size-convergence fit parameters, curves, errors, and flags.');

score = round(score, 6);
report = struct( ...
    'record_id', 'matlab_gpumd_conductivity_robustness', ...
    'score', score, ...
    'max_score', maxScore, ...
    'passed', score >= maxScore - 1.0e-9, ...
    'details', details);

reportDirectory = fileparts(report_path);
if strlength(reportDirectory) > 0 && ~isfolder(reportDirectory)
    mkdir(reportDirectory);
end
jsonText = jsonencode(report, 'PrettyPrint', true);
fileId = fopen(report_path, 'w', 'n', 'UTF-8');
if fileId < 0
    error('GPUMD:WriteReport', 'Cannot write evaluation report: %s', report_path);
end
fileCleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, '%s\n', jsonText);
clear fileCleanup;
end

function tf = hasColumns(inputTable, requiredNames)
tf = all(ismember(requiredNames, string(inputTable.Properties.VariableNames)));
end

function matches = numericMatches(candidate, truth, tolerance)
try
    candidate = numericVector(candidate);
    truth = numericVector(truth);
    matches = isfinite(candidate) & isfinite(truth) & ...
        abs(candidate - truth) <= tolerance;
catch
    matches = false(size(truth));
end
end

function values = numericVector(inputValues)
if isnumeric(inputValues) || islogical(inputValues)
    values = double(inputValues);
else
    values = str2double(string(inputValues));
end
end

function tf = scalarMatch(candidate, truth, tolerance)
try
    candidate = double(candidate);
    truth = double(truth);
    tf = isscalar(candidate) && isscalar(truth) && isfinite(candidate) && ...
        isfinite(truth) && abs(candidate - truth) <= tolerance;
catch
    tf = false;
end
end

function tf = arrayMatch(candidate, truth, tolerance)
try
    candidate = numericVector(candidate);
    truth = numericVector(truth);
    tf = numel(candidate) == numel(truth) && all(isfinite(candidate), 'all') && ...
        all(isfinite(truth), 'all') && ...
        all(abs(candidate(:) - truth(:)) <= tolerance);
catch
    tf = false;
end
end

function earned = robustnessScore(candidateSummary, truthSummary, sectionName, ...
        materialNames, groupLabels)
earned = 0.0;
if ~isstruct(candidateSummary) || ~isfield(candidateSummary, sectionName) || ...
        ~isfield(truthSummary, sectionName)
    return;
end

candidateSection = candidateSummary.(sectionName);
truthSection = truthSummary.(sectionName);
if ~isstruct(candidateSection)
    return;
end
pointsPerMaterial = 10.0 / numel(materialNames);

for i = 1:numel(materialNames)
    material = char(materialNames(i));
    if ~isfield(candidateSection, material) || ~isfield(truthSection, material)
        continue;
    end
    candidateEntry = candidateSection.(material);
    truthEntry = truthSection.(material);
    if ~isstruct(candidateEntry)
        continue;
    end
    correct = 0;
    checks = numel(groupLabels) + 2;

    if isfield(candidateEntry, 'group_means_W_mK') && ...
            isstruct(candidateEntry.group_means_W_mK) && ...
            isfield(truthEntry, 'group_means_W_mK')
        for j = 1:numel(groupLabels)
            label = char(groupLabels(j));
            if isfield(candidateEntry.group_means_W_mK, label) && ...
                    isfield(truthEntry.group_means_W_mK, label) && ...
                    scalarMatch(candidateEntry.group_means_W_mK.(label), ...
                    truthEntry.group_means_W_mK.(label), 1.0e-6)
                correct = correct + 1;
            end
        end
    end

    if isfield(candidateEntry, 'max_relative_deviation') && ...
            isfield(truthEntry, 'max_relative_deviation') && ...
            scalarMatch(candidateEntry.max_relative_deviation, ...
            truthEntry.max_relative_deviation, 1.0e-6)
        correct = correct + 1;
    end
    if isfield(candidateEntry, 'robust_within_5_percent') && ...
            isfield(truthEntry, 'robust_within_5_percent') && ...
            isequal(candidateEntry.robust_within_5_percent, ...
            truthEntry.robust_within_5_percent)
        correct = correct + 1;
    end

    earned = earned + pointsPerMaterial * correct / checks;
end
end

function earned = sizeConvergenceScore(candidateSummary, truthSummary, materialNames)
earned = 0.0;
sectionName = 'size_convergence_fit';
if ~isstruct(candidateSummary) || ~isfield(candidateSummary, sectionName) || ...
        ~isfield(truthSummary, sectionName)
    return;
end

candidateSection = candidateSummary.(sectionName);
truthSection = truthSummary.(sectionName);
if ~isstruct(candidateSection)
    return;
end

pointsPerMaterial = 12.0 / numel(materialNames);
numericFields = ["k_inf_W_mK", "a_W_nm_mK", "rmse_W_mK", "relative_rmse"];
tolerances = [1.0e-6, 1.0e-5, 1.0e-7, 1.0e-7];
checksPerMaterial = numel(numericFields) + 3;

for i = 1:numel(materialNames)
    material = char(materialNames(i));
    if ~isfield(candidateSection, material) || ~isfield(truthSection, material)
        continue;
    end
    candidateEntry = candidateSection.(material);
    truthEntry = truthSection.(material);
    if ~isstruct(candidateEntry)
        continue;
    end

    correct = 0;
    for j = 1:numel(numericFields)
        field = char(numericFields(j));
        if isfield(candidateEntry, field) && isfield(truthEntry, field) && ...
                scalarMatch(candidateEntry.(field), truthEntry.(field), tolerances(j))
            correct = correct + 1;
        end
    end
    if isfield(candidateEntry, 'observed_group_means_W_mK') && ...
            isfield(truthEntry, 'observed_group_means_W_mK') && ...
            arrayMatch(candidateEntry.observed_group_means_W_mK, ...
            truthEntry.observed_group_means_W_mK, 1.0e-6)
        correct = correct + 1;
    end
    if isfield(candidateEntry, 'fitted_group_means_W_mK') && ...
            isfield(truthEntry, 'fitted_group_means_W_mK') && ...
            arrayMatch(candidateEntry.fitted_group_means_W_mK, ...
            truthEntry.fitted_group_means_W_mK, 1.0e-6)
        correct = correct + 1;
    end
    if isfield(candidateEntry, 'converged') && isfield(truthEntry, 'converged') && ...
            isequal(candidateEntry.converged, truthEntry.converged)
        correct = correct + 1;
    end

    earned = earned + pointsPerMaterial * correct / checksPerMaterial;
end
end

function [score, details] = award(score, details, id, points, earned, message)
earned = max(0.0, min(points, earned));
score = score + earned;
details(end + 1) = struct( ...
    'id', id, ...
    'points', points, ...
    'earned', round(earned, 6), ...
    'message', message);
end
