function report = evaluate_submission(submission_dir, ground_truth_dir, report_path)
%EVALUATE_SUBMISSION Deterministically grade VACF/PDOS signal outputs.
%   REPORT = EVALUATE_SUBMISSION(SUBMISSION_DIR, GROUND_TRUTH_DIR, REPORT_PATH)
%   returns a report struct and writes the same report as UTF-8 JSON.

arguments
    submission_dir (1, 1) string
    ground_truth_dir (1, 1) string
    report_path (1, 1) string
end

materialNames = ["CsSnBr3", "Cs2SnBr6", "hetero_N2"];
maxScore = 100.0;
score = 0.0;
details = struct('id', {}, 'points', {}, 'earned', {}, 'message', {});

truthVacfPath = fullfile(ground_truth_dir, 'vacf.csv');
truthPdosPath = fullfile(ground_truth_dir, 'pdos.csv');
truthSummaryPath = fullfile(ground_truth_dir, 'phonon_summary.json');
if ~isfile(truthVacfPath) || ~isfile(truthPdosPath) || ~isfile(truthSummaryPath)
    error('VACFPDOS:MissingGroundTruth', ...
        'Ground Truth files are missing from %s.', ground_truth_dir);
end

truthVacf = readtable(truthVacfPath, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
truthPdos = readtable(truthPdosPath, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
truthSummary = jsondecode(fileread(truthSummaryPath));

candidateVacfPath = fullfile(submission_dir, 'vacf.csv');
candidatePdosPath = fullfile(submission_dir, 'pdos.csv');
candidateSummaryPath = fullfile(submission_dir, 'phonon_summary.json');
vacfExists = isfile(candidateVacfPath);
pdosExists = isfile(candidatePdosPath);
summaryExists = isfile(candidateSummaryPath);

vacfParsed = false;
pdosParsed = false;
summaryParsed = false;
candidateVacf = table();
candidatePdos = table();
candidateSummary = struct();

if vacfExists
    try
        candidateVacf = readtable(candidateVacfPath, 'FileType', 'text', ...
            'Delimiter', ',', 'TextType', 'string');
        vacfParsed = true;
    catch
        vacfParsed = false;
    end
end
if pdosExists
    try
        candidatePdos = readtable(candidatePdosPath, 'FileType', 'text', ...
            'Delimiter', ',', 'TextType', 'string');
        pdosParsed = true;
    catch
        pdosParsed = false;
    end
end
if summaryExists
    try
        candidateSummary = jsondecode(fileread(candidateSummaryPath));
        summaryParsed = true;
    catch
        summaryParsed = false;
    end
end

earned = double(vacfExists) + double(pdosExists) + double(summaryExists) + ...
    double(vacfParsed) + double(pdosParsed) + double(summaryParsed);
if summaryParsed && isfield(candidateSummary, 'schema_version') && ...
        strcmp(string(candidateSummary.schema_version), "1.0")
    earned = earned + 1.0;
end
if summaryParsed && isfield(candidateSummary, 'synthetic_data') && ...
        isequal(candidateSummary.synthetic_data, true)
    earned = earned + 1.0;
end
if summaryParsed && isfield(candidateSummary, 'systems') && ...
        all(isfield(candidateSummary.systems, cellstr(materialNames)))
    earned = earned + 1.0;
end
requiredTop = {'interface_statistics', 'comparisons', 'physical_conclusion', 'method'};
if summaryParsed && all(isfield(candidateSummary, requiredTop))
    earned = earned + 1.0;
end
[score, details] = award(score, details, 'files_and_schema', 10.0, earned, ...
    'Required files, parsing, schema declaration, and system keys.');

requiredVacf = ["material", "time_ps", "vacf"];
requiredPdos = ["material", "frequency_THz", "pdos"];
vacfReady = vacfParsed && hasColumns(candidateVacf, requiredVacf) && ...
    hasColumns(truthVacf, requiredVacf);
pdosReady = pdosParsed && hasColumns(candidatePdos, requiredPdos) && ...
    hasColumns(truthPdos, requiredPdos);

vacfCoverage = zeros(1, numel(materialNames));
vacfCorrelation = -ones(1, numel(materialNames));
pdosCoverage = zeros(1, numel(materialNames));
pdosCorrelation = -ones(1, numel(materialNames));
pdosNonnegative = false(1, numel(materialNames));
for index = 1:numel(materialNames)
    material = materialNames(index);
    if vacfReady
        [vacfCoverage(index), vacfCorrelation(index)] = compareCurve( ...
            candidateVacf, truthVacf, material, 'time_ps', 'vacf');
    end
    if pdosReady
        [pdosCoverage(index), pdosCorrelation(index), pdosNonnegative(index)] = ...
            compareCurve(candidatePdos, truthPdos, material, ...
            'frequency_THz', 'pdos');
    end
end

earned = 5.0 * mean(vacfCoverage) + 5.0 * mean(pdosCoverage);
[score, details] = award(score, details, 'curve_coverage', 10.0, earned, ...
    'Rounded-axis coverage for all materials and both curve files.');

earned = 0.0;
for index = 1:numel(materialNames)
    if vacfCoverage(index) >= 0.95 && vacfCorrelation(index) >= 0.995
        earned = earned + 5.0;
    elseif vacfCoverage(index) >= 0.80 && vacfCorrelation(index) >= 0.98
        earned = earned + 2.5;
    end
end
[score, details] = award(score, details, 'vacf_curves', 15.0, earned, ...
    'Unbiased normalized VACF curve coverage and correlation.');

earned = 0.0;
for index = 1:numel(materialNames)
    if pdosNonnegative(index) && pdosCoverage(index) >= 0.95 && ...
            pdosCorrelation(index) >= 0.995
        earned = earned + 5.0;
    elseif pdosNonnegative(index) && pdosCoverage(index) >= 0.80 && ...
            pdosCorrelation(index) >= 0.98
        earned = earned + 2.5;
    end
end
[score, details] = award(score, details, 'pdos_curves', 15.0, earned, ...
    'Symmetric-Hann PDOS curve coverage, correlation, and non-negativity.');

earned = 0.0;
if vacfReady
    for index = 1:numel(materialNames)
        rows = string(candidateVacf.material) == materialNames(index);
        time = numericVector(candidateVacf.time_ps(rows));
        values = numericVector(candidateVacf.vacf(rows));
        zeroLocation = find(abs(time) <= 1.0e-8, 1, 'first');
        if ~isempty(zeroLocation) && isfinite(values(zeroLocation)) && ...
                abs(values(zeroLocation) - 1.0) <= 1.0e-6
            earned = earned + 1.0;
        end
    end
end
if pdosReady
    for index = 1:numel(materialNames)
        rows = string(candidatePdos.material) == materialNames(index);
        frequency = numericVector(candidatePdos.frequency_THz(rows));
        values = numericVector(candidatePdos.pdos(rows));
        [frequency, order] = sort(frequency);
        values = values(order);
        if ~isempty(values) && all(isfinite(values)) && all(values >= -1.0e-10)
            earned = earned + 1.0;
        end
        if numel(frequency) >= 2 && all(isfinite(frequency)) && ...
                all(diff(frequency) > 0) && ...
                abs(trapz(frequency, values) - 1.0) <= 1.0e-3
            earned = earned + 1.0;
        end
    end
end
if vacfReady && pdosReady && axesValid(candidateVacf, 'time_ps', materialNames) && ...
        axesValid(candidatePdos, 'frequency_THz', materialNames)
    earned = earned + 1.0;
end
[score, details] = award(score, details, 'normalization_and_axes', 10.0, earned, ...
    'VACF lag-zero normalization, PDOS invariants, and valid axes.');

earned = 0.0;
if summaryParsed && isfield(candidateSummary, 'systems') && ...
        isfield(truthSummary, 'systems')
    for index = 1:numel(materialNames)
        material = char(materialNames(index));
        if summaryScalarMatch(candidateSummary, truthSummary, material, ...
                'vacf_positive_integral_ps', 0.02)
            earned = earned + 1.5;
        end
        if summaryScalarMatch(candidateSummary, truthSummary, material, ...
                'band_fraction_2_3_THz', 0.02)
            earned = earned + 1.5;
        end
        if summaryScalarMatch(candidateSummary, truthSummary, material, ...
                'dominant_peak_THz', 0.11)
            earned = earned + 0.5;
        end
        if summaryScalarMatch(candidateSummary, truthSummary, material, ...
                'band_peak_THz', 0.11)
            earned = earned + 0.5;
        end
    end
    if summaryScalarMatch(candidateSummary, truthSummary, 'hetero_N2', ...
            'interface_share_2_3_THz', 0.03)
        earned = earned + 3.0;
    end
end
[score, details] = award(score, details, 'summary_metrics', 15.0, earned, ...
    'VACF integrals, 2-3 THz fractions, interface share, and peaks.');

earned = 0.0;
if summaryParsed && isfield(candidateSummary, 'interface_statistics') && ...
        isfield(truthSummary, 'interface_statistics')
    candidateStats = candidateSummary.interface_statistics;
    truthStats = truthSummary.interface_statistics;
    if isfield(candidateStats, 'test') && isfield(truthStats, 'test') && ...
            strcmp(string(candidateStats.test), string(truthStats.test))
        earned = earned + 1.0;
    end
    if statisticsScalarMatch(candidateStats, truthStats, 'alpha', 1.0e-12, 0.0)
        earned = earned + 1.0;
    end
    for field = ["n_interface", "n_noninterface"]
        if statisticsScalarMatch(candidateStats, truthStats, char(field), 0.0, 0.0)
            earned = earned + 0.5;
        end
    end
    for field = ["mean_interface_band_power", "mean_noninterface_band_power", ...
            "mean_difference", "ci95_low", "ci95_high"]
        if statisticsScalarMatch(candidateStats, truthStats, char(field), 1.0e-12, 1.0e-6)
            earned = earned + 1.0;
        end
    end
    for field = ["t_statistic", "degrees_of_freedom", "p_value"]
        if statisticsScalarMatch(candidateStats, truthStats, char(field), 1.0e-12, 1.0e-5)
            earned = earned + 1.0;
        end
    end
    if isfield(candidateStats, 'reject_equal_means') && ...
            isfield(truthStats, 'reject_equal_means') && ...
            isequal(candidateStats.reject_equal_means, truthStats.reject_equal_means)
        earned = earned + 2.0;
    end
    if statisticsScalarMatch(candidateStats, truthStats, 'hedges_g', 1.0e-12, 1.0e-5)
        earned = earned + 2.0;
    end
end
[score, details] = award(score, details, 'interface_statistics', 15.0, earned, ...
    'Welch test, confidence interval, and Hedges g for spatial channel groups.');

earned = 0.0;
comparisonFields = ["heterostructure_shorter_vacf_than_bulks", ...
    "heterostructure_enhanced_2_3_THz_fraction", ...
    "heterostructure_interface_localized_2_3_THz", ...
    "heterostructure_interface_band_power_significant"];
if summaryParsed && isfield(candidateSummary, 'comparisons') && ...
        isfield(truthSummary, 'comparisons')
    for index = 1:numel(comparisonFields)
        field = char(comparisonFields(index));
        if isfield(candidateSummary.comparisons, field) && ...
                isfield(truthSummary.comparisons, field) && ...
                isequal(candidateSummary.comparisons.(field), ...
                truthSummary.comparisons.(field))
            earned = earned + 2.0;
        end
    end
end
if summaryParsed && isfield(candidateSummary, 'physical_conclusion') && ...
        isfield(truthSummary, 'physical_conclusion') && ...
        strcmp(string(candidateSummary.physical_conclusion), ...
        string(truthSummary.physical_conclusion))
    earned = earned + 1.0;
end
if summaryParsed && isfield(candidateSummary, 'method') && ...
        isfield(truthSummary, 'method') && ...
        strcmp(string(candidateSummary.method), string(truthSummary.method))
    earned = earned + 1.0;
end
[score, details] = award(score, details, 'physical_conclusions', 10.0, earned, ...
    'Four threshold decisions, fixed conclusion, and declared method.');

score = round(score, 6);
report = struct( ...
    'record_id', 'matlab_vacf_pdos_interface_localization', ...
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
    error('VACFPDOS:WriteReport', ...
        'Cannot write evaluation report: %s', report_path);
end
fileCleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, '%s\n', jsonText);
clear fileCleanup;
end

function tf = hasColumns(inputTable, requiredNames)
tf = all(ismember(requiredNames, string(inputTable.Properties.VariableNames)));
end

function [coverage, correlationValue, nonnegative] = compareCurve( ...
        candidate, truth, material, xName, yName)
candidateRows = string(candidate.material) == material;
truthRows = string(truth.material) == material;
candidateX = numericVector(candidate.(xName)(candidateRows));
candidateY = numericVector(candidate.(yName)(candidateRows));
truthX = numericVector(truth.(xName)(truthRows));
truthY = numericVector(truth.(yName)(truthRows));
nonnegative = ~isempty(candidateY) && all(isfinite(candidateY)) && ...
    all(candidateY >= -1.0e-10);

if isempty(truthX) || isempty(candidateX)
    coverage = 0.0;
    correlationValue = -1.0;
    return;
end

candidateKeys = round(candidateX * 1.0e8) / 1.0e8;
truthKeys = round(truthX * 1.0e8) / 1.0e8;
[candidateKeys, uniqueLocations] = unique(candidateKeys, 'stable');
candidateY = candidateY(uniqueLocations);
[matched, candidateLocations] = ismember(truthKeys, candidateKeys);
valid = matched & isfinite(truthY);
matchedLocations = candidateLocations(valid);
valid(valid) = isfinite(candidateY(matchedLocations));
coverage = nnz(valid) / numel(truthKeys);
correlationValue = safeCorrelation(candidateY(candidateLocations(valid)), truthY(valid));
end

function value = safeCorrelation(left, right)
left = left(:);
right = right(:);
if numel(left) ~= numel(right) || numel(left) < 2 || ...
        any(~isfinite(left)) || any(~isfinite(right))
    value = -1.0;
    return;
end
left = left - mean(left);
right = right - mean(right);
denominator = norm(left) * norm(right);
if denominator <= 0
    value = -1.0;
else
    value = (left' * right) / denominator;
end
end

function tf = axesValid(inputTable, xName, materialNames)
tf = true;
for index = 1:numel(materialNames)
    rows = string(inputTable.material) == materialNames(index);
    values = numericVector(inputTable.(xName)(rows));
    if numel(values) < 2 || any(~isfinite(values)) || ...
            numel(unique(round(values * 1.0e8) / 1.0e8)) ~= numel(values)
        tf = false;
        return;
    end
end
end

function tf = summaryScalarMatch(candidate, truth, material, field, tolerance)
tf = false;
if ~isfield(candidate.systems, material) || ~isfield(truth.systems, material)
    return;
end
candidateEntry = candidate.systems.(material);
truthEntry = truth.systems.(material);
if ~isfield(candidateEntry, field) || ~isfield(truthEntry, field)
    return;
end
try
    candidateValue = double(candidateEntry.(field));
    truthValue = double(truthEntry.(field));
    tf = isscalar(candidateValue) && isscalar(truthValue) && ...
        isfinite(candidateValue) && isfinite(truthValue) && ...
        abs(candidateValue - truthValue) <= tolerance;
catch
    tf = false;
end
end

function tf = statisticsScalarMatch(candidate, truth, field, absoluteTolerance, ...
        relativeTolerance)
tf = false;
if ~isfield(candidate, field) || ~isfield(truth, field)
    return;
end
try
    candidateValue = double(candidate.(field));
    truthValue = double(truth.(field));
    tolerance = absoluteTolerance + relativeTolerance * abs(truthValue);
    tf = isscalar(candidateValue) && isscalar(truthValue) && ...
        isfinite(candidateValue) && isfinite(truthValue) && ...
        abs(candidateValue - truthValue) <= tolerance;
catch
    tf = false;
end
end

function values = numericVector(values)
if isnumeric(values) || islogical(values)
    values = double(values(:));
else
    values = str2double(string(values(:)));
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
