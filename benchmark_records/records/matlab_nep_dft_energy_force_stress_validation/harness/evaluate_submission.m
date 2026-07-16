function report = evaluate_submission(submission_dir, ground_truth_dir, report_path)
%EVALUATE_SUBMISSION Deterministically score one benchmark submission.

submission_dir = char(submission_dir);
ground_truth_dir = char(ground_truth_dir);
report_path = char(report_path);
details = struct('criterion', {}, 'max_points', {}, ...
    'earned_points', {}, 'passed_checks', {}, 'total_checks', {}, ...
    'message', {});
score = 0.0;

required_files = {'metrics.json', 'residuals.csv', ...
    'training_assessment.json', 'calibration.json'};
submission_paths = cellfun(@(name) fullfile(submission_dir, name), ...
    required_files, 'UniformOutput', false);
ground_truth_paths = cellfun(@(name) fullfile(ground_truth_dir, name), ...
    required_files, 'UniformOutput', false);
submission_exists = cellfun(@isfile, submission_paths);
ground_truth_exists = cellfun(@isfile, ground_truth_paths);

if ~all(ground_truth_exists)
    error('Ground Truth directory is incomplete: %s', ground_truth_dir);
end

if ~all(submission_exists)
    passed = sum(submission_exists);
    earned = 8.0 * passed / numel(required_files);
    details = appendDetail(details, 'required_files_and_schema', ...
        8.0, earned, passed, numel(required_files), ...
        'One or more required submission files are missing.');
    score = score + earned;
    report = finalizeReport(score, details, report_path);
    return;
end

try
    actual_metrics = jsondecode(fileread(submission_paths{1}));
    actual_residuals = readtable(submission_paths{2}, ...
        'TextType', 'string', 'VariableNamingRule', 'preserve');
    actual_training = jsondecode(fileread(submission_paths{3}));
    actual_calibration = jsondecode(fileread(submission_paths{4}));
    expected_metrics = jsondecode(fileread(ground_truth_paths{1}));
    expected_residuals = readtable(ground_truth_paths{2}, ...
        'TextType', 'string', 'VariableNamingRule', 'preserve');
    expected_training = jsondecode(fileread(ground_truth_paths{3}));
    expected_calibration = jsondecode(fileread(ground_truth_paths{4}));
catch exception
    details = appendDetail(details, 'required_files_and_schema', ...
        8.0, 0.0, 0, 1, ['Parse failure: ' exception.message]);
    report = finalizeReport(score, details, report_path);
    return;
end

schema_checks = false(10, 1);
schema_checks(1) = textFieldEquals(actual_metrics, ...
    'schema_version', '1.0');
schema_checks(2) = trueField(actual_metrics, 'synthetic_data');
schema_checks(3) = textFieldEquals(actual_training, ...
    'schema_version', '1.0');
schema_checks(4) = trueField(actual_training, 'synthetic_data');
schema_checks(5) = isfield(actual_metrics, 'overall');
schema_checks(6) = all(isfield(actual_metrics, ...
    {'CsSnBr3', 'Cs2SnBr6', 'hetero_N2'}));
schema_checks(7) = stringListsEqual(getFieldOrEmpty(actual_metrics, ...
    'stress_voigt_order'), {'xx', 'yy', 'zz', 'yz', 'xz', 'xy'});
schema_checks(8) = textFieldEquals(actual_calibration, ...
    'schema_version', '1.0');
schema_checks(9) = trueField(actual_calibration, 'synthetic_data');
schema_checks(10) = textFieldEquals(actual_calibration, ...
    'method', 'bounded_regularized_lsqlin');
[details, earned] = scoreChecks(details, ...
    'required_files_and_schema', 8.0, schema_checks, ...
    'Required files parsed and schema markers were checked.');
score = score + earned;

[passed, total] = compareMetricFields(actual_metrics, expected_metrics, ...
    'overall', {'count', 'energy_RMSE_meV_atom', 'energy_R2'}, ...
    [0, 1e-6, 1e-8], 1e-6);
[details, earned] = scoreCount(details, 'overall_energy_metrics', ...
    8.0, passed, total, ...
    'Overall count, energy RMSE, and energy R2 were compared.');
score = score + earned;

[passed, total] = compareMetricFields(actual_metrics, expected_metrics, ...
    'overall', {'force_RMSE_eV_A', 'force_R2'}, ...
    [1e-8, 1e-8], 1e-6);
[details, earned] = scoreCount(details, 'overall_force_metrics', ...
    12.0, passed, total, ...
    'Overall force RMSE and force R2 were compared.');
score = score + earned;

[passed, total] = compareMetricFields(actual_metrics, expected_metrics, ...
    'overall', {'stress_RMSE_GPa', 'stress_R2'}, ...
    [1e-8, 1e-8], 1e-6);
[details, earned] = scoreCount(details, 'overall_stress_metrics', ...
    8.0, passed, total, ...
    'Overall stress RMSE and stress R2 were compared.');
score = score + earned;

material_groups = {'CsSnBr3', 'Cs2SnBr6', 'hetero_N2'};
metric_fields = {'count', 'energy_RMSE_meV_atom', 'energy_R2', ...
    'force_RMSE_eV_A', 'force_R2', 'stress_RMSE_GPa', 'stress_R2'};
metric_tolerances = [0, 1e-6, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8];
passed = 0;
total = 0;
for group_index = 1:numel(material_groups)
    [group_passed, group_total] = compareMetricFields( ...
        actual_metrics, expected_metrics, material_groups{group_index}, ...
        metric_fields, metric_tolerances, 1e-6);
    passed = passed + group_passed;
    total = total + group_total;
end
[details, earned] = scoreCount(details, 'per_material_metrics', ...
    16.0, passed, total, ...
    'All six statistics and count were checked for each material group.');
score = score + earned;

residual_columns = {'structure_id', 'material', 'n_atoms', ...
    'energy_ref_eV', 'energy_pred_eV', 'energy_error_meV_atom', ...
    'force_rmse_eV_A', 'stress_rmse_GPa', ...
    'max_force_abs_error_eV_A'};
columns_ok = all(ismember(residual_columns, ...
    actual_residuals.Properties.VariableNames));
expected_columns_ok = all(ismember(residual_columns, ...
    expected_residuals.Properties.VariableNames));

coverage_checks = false(5, 1);
aligned_actual = table;
if columns_ok && expected_columns_ok
    actual_ids = string(actual_residuals.structure_id);
    expected_ids = string(expected_residuals.structure_id);
    [matched, locations] = ismember(expected_ids, actual_ids);
    coverage_checks(1) = height(actual_residuals) == ...
        height(expected_residuals);
    coverage_checks(2) = numel(unique(actual_ids)) == numel(actual_ids);
    coverage_checks(3) = all(matched) && numel(actual_ids) == ...
        numel(expected_ids);
    if all(matched)
        aligned_actual = actual_residuals(locations, :);
        coverage_checks(4) = all(string(aligned_actual.material) == ...
            string(expected_residuals.material));
        coverage_checks(5) = all(double(aligned_actual.n_atoms) == ...
            double(expected_residuals.n_atoms));
    end
end
[details, earned] = scoreChecks(details, ...
    'structure_and_residual_coverage', 8.0, coverage_checks, ...
    'Residual row count, unique IDs, exact ID set, material, and atom count were checked.');
score = score + earned;

residual_passed = 0;
residual_total = 1;
if ~isempty(aligned_actual)
    residual_numeric_fields = {'energy_ref_eV', 'energy_pred_eV', ...
        'energy_error_meV_atom', 'force_rmse_eV_A', ...
        'stress_rmse_GPa', 'max_force_abs_error_eV_A'};
    residual_passed = 0;
    residual_total = 0;
    for field_index = 1:numel(residual_numeric_fields)
        field_name = residual_numeric_fields{field_index};
        actual_values = double(aligned_actual.(field_name));
        expected_values = double(expected_residuals.(field_name));
        for value_index = 1:numel(expected_values)
            residual_total = residual_total + 1;
            if numericClose(actual_values(value_index), ...
                    expected_values(value_index), 1e-7, 1e-6)
                residual_passed = residual_passed + 1;
            end
        end
    end
end
[details, earned] = scoreCount(details, ...
    'per_structure_residual_values', 12.0, residual_passed, ...
    residual_total, ...
    'Per-structure energy, force, and stress residual values were compared by structure_id.');
score = score + earned;

training_passed = 0;
training_total = 0;
exact_training_fields = {'n_records', 'first_step', 'last_step', ...
    'best_test_step'};
numeric_training_fields = {'initial_train_total_loss', ...
    'final_train_total_loss', 'initial_test_total_loss', ...
    'final_test_total_loss', 'train_loss_reduction_fraction', ...
    'test_loss_reduction_fraction', ...
    'final_to_best_test_loss_ratio', ...
    'final_test_to_train_loss_ratio'};
boolean_training_fields = {'converged', 'no_overfitting', ...
    'metric_quality_pass', 'reliable_for_thermal_transport'};

for field_index = 1:numel(exact_training_fields)
    field_name = exact_training_fields{field_index};
    training_total = training_total + 1;
    if numericFieldsEqual(actual_training, expected_training, ...
            field_name, 0, 0)
        training_passed = training_passed + 1;
    end
end
for field_index = 1:numel(numeric_training_fields)
    field_name = numeric_training_fields{field_index};
    training_total = training_total + 1;
    if numericFieldsEqual(actual_training, expected_training, ...
            field_name, 1e-9, 1e-6)
        training_passed = training_passed + 1;
    end
end
for field_index = 1:numel(boolean_training_fields)
    field_name = boolean_training_fields{field_index};
    training_total = training_total + 1;
    if booleanFieldsEqual(actual_training, expected_training, field_name)
        training_passed = training_passed + 1;
    end
end
[details, earned] = scoreCount(details, 'training_and_reliability', ...
    10.0, training_passed, training_total, ...
    'Training/test loss statistics and reliability decisions were compared.');
score = score + earned;

calibration_passed = 0;
calibration_total = 0;
calibration_passed = calibration_passed + textFieldEquals( ...
    actual_calibration, 'method', 'bounded_regularized_lsqlin');
calibration_total = calibration_total + 1;

if isfield(actual_calibration, 'objective') && ...
        isfield(expected_calibration, 'objective')
    objective_fields = {'energy_scale_eV_atom', 'force_scale_eV_A', ...
        'stress_scale_GPa', 'regularization_lambda'};
    for field_index = 1:numel(objective_fields)
        field_name = objective_fields{field_index};
        calibration_total = calibration_total + 1;
        if numericFieldsEqual(actual_calibration.objective, ...
                expected_calibration.objective, field_name, 1e-12, 1e-9)
            calibration_passed = calibration_passed + 1;
        end
    end
    calibration_total = calibration_total + 1;
    if booleanFieldsEqual(actual_calibration.objective, ...
            expected_calibration.objective, 'equal_block_weighting')
        calibration_passed = calibration_passed + 1;
    end
    objective_array_fields = {'regularization_center', ...
        'regularization_scales'};
    for field_index = 1:numel(objective_array_fields)
        field_name = objective_array_fields{field_index};
        calibration_total = calibration_total + 1;
        if numericArrayFieldsEqual(actual_calibration.objective, ...
                expected_calibration.objective, field_name, 1e-12, 1e-9)
            calibration_passed = calibration_passed + 1;
        end
    end
else
    calibration_total = calibration_total + 7;
end

if isfield(actual_calibration, 'bounds') && ...
        isfield(expected_calibration, 'bounds')
    bound_fields = {'energy_bias_eV_atom', 'force_scale', 'stress_scale'};
    for field_index = 1:numel(bound_fields)
        field_name = bound_fields{field_index};
        calibration_total = calibration_total + 1;
        if numericArrayFieldsEqual(actual_calibration.bounds, ...
                expected_calibration.bounds, field_name, 1e-12, 1e-9)
            calibration_passed = calibration_passed + 1;
        end
    end
else
    calibration_total = calibration_total + 3;
end

if isfield(actual_calibration, 'parameters') && ...
        isfield(expected_calibration, 'parameters')
    parameter_fields = {'energy_bias_eV_atom', 'force_scale', ...
        'stress_scale'};
    for field_index = 1:numel(parameter_fields)
        field_name = parameter_fields{field_index};
        calibration_total = calibration_total + 1;
        if numericFieldsEqual(actual_calibration.parameters, ...
                expected_calibration.parameters, field_name, 1e-9, 1e-7)
            calibration_passed = calibration_passed + 1;
        end
    end
else
    calibration_total = calibration_total + 3;
end

if isfield(actual_calibration, 'solver') && ...
        isfield(expected_calibration, 'solver')
    calibration_total = calibration_total + 1;
    if textFieldEquals(actual_calibration.solver, 'name', 'lsqlin')
        calibration_passed = calibration_passed + 1;
    end
    solver_fields = {'exitflag', 'iterations', 'resnorm'};
    solver_tolerances = [0, 0, 1e-9];
    for field_index = 1:numel(solver_fields)
        field_name = solver_fields{field_index};
        calibration_total = calibration_total + 1;
        if numericFieldsEqual(actual_calibration.solver, ...
                expected_calibration.solver, field_name, ...
                solver_tolerances(field_index), 1e-7)
            calibration_passed = calibration_passed + 1;
        end
    end
else
    calibration_total = calibration_total + 4;
end

calibration_total = calibration_total + 1;
if booleanFieldsEqual(actual_calibration, expected_calibration, ...
        'unique_solution')
    calibration_passed = calibration_passed + 1;
end
if isfield(actual_calibration, 'strict_convexity') && ...
        isfield(expected_calibration, 'strict_convexity')
    convex_fields = {'augmented_design_rank', 'parameter_count', ...
        'minimum_gram_eigenvalue'};
    convex_tolerances = [0, 0, 1e-9];
    for field_index = 1:numel(convex_fields)
        field_name = convex_fields{field_index};
        calibration_total = calibration_total + 1;
        if numericFieldsEqual(actual_calibration.strict_convexity, ...
                expected_calibration.strict_convexity, field_name, ...
                convex_tolerances(field_index), 1e-7)
            calibration_passed = calibration_passed + 1;
        end
    end
    calibration_total = calibration_total + 1;
    if booleanFieldsEqual(actual_calibration.strict_convexity, ...
            expected_calibration.strict_convexity, 'verified')
        calibration_passed = calibration_passed + 1;
    end
else
    calibration_total = calibration_total + 4;
end
[details, earned] = scoreCount(details, ...
    'constrained_calibration_parameters', 8.0, calibration_passed, ...
    calibration_total, ...
    'Fixed objective, bounds, regularization, lsqlin solution, and strict convexity were checked.');
score = score + earned;

calibrated_metric_passed = 0;
calibrated_metric_total = 0;
phase_names = {'metrics_before', 'metrics_after'};
all_metric_fields = {'count', 'energy_RMSE_meV_atom', 'energy_R2', ...
    'force_RMSE_eV_A', 'force_R2', 'stress_RMSE_GPa', 'stress_R2'};
all_metric_tolerances = [0, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8];
all_groups = {'overall', 'CsSnBr3', 'Cs2SnBr6', 'hetero_N2'};
for phase_index = 1:numel(phase_names)
    phase_name = phase_names{phase_index};
    if isfield(actual_calibration, phase_name) && ...
            isfield(expected_calibration, phase_name)
        actual_phase = actual_calibration.(phase_name);
        expected_phase = expected_calibration.(phase_name);
        for group_index = 1:numel(all_groups)
            [group_passed, group_total] = compareMetricFields( ...
                actual_phase, expected_phase, all_groups{group_index}, ...
                all_metric_fields, all_metric_tolerances, 1e-6);
            calibrated_metric_passed = calibrated_metric_passed + ...
                group_passed;
            calibrated_metric_total = calibrated_metric_total + group_total;
        end
    else
        calibrated_metric_total = calibrated_metric_total + ...
            numel(all_groups) * numel(all_metric_fields);
    end
end
[details, earned] = scoreCount(details, ...
    'calibrated_before_after_metrics', 8.0, calibrated_metric_passed, ...
    calibrated_metric_total, ...
    'Before/after grouped RMSE and R2 values were compared.');
score = score + earned;

distribution_passed = 0;
distribution_total = 0;
distribution_phases = {'residual_statistics_before', ...
    'residual_statistics_after'};
distribution_groups = {'energy_meV_atom', 'force_eV_A', 'stress_GPa'};
distribution_fields = {'mean', 'std', 'rmse', 'count'};
distribution_tolerances = [1e-8, 1e-8, 1e-8, 0];
for phase_index = 1:numel(distribution_phases)
    phase_name = distribution_phases{phase_index};
    if isfield(actual_calibration, phase_name) && ...
            isfield(expected_calibration, phase_name)
        actual_phase = actual_calibration.(phase_name);
        expected_phase = expected_calibration.(phase_name);
        for group_index = 1:numel(distribution_groups)
            group_name = distribution_groups{group_index};
            if isfield(actual_phase, group_name) && ...
                    isfield(expected_phase, group_name)
                for field_index = 1:numel(distribution_fields)
                    distribution_total = distribution_total + 1;
                    if numericFieldsEqual(actual_phase.(group_name), ...
                            expected_phase.(group_name), ...
                            distribution_fields{field_index}, ...
                            distribution_tolerances(field_index), 1e-6)
                        distribution_passed = distribution_passed + 1;
                    end
                end
            else
                distribution_total = distribution_total + ...
                    numel(distribution_fields);
            end
        end
    else
        distribution_total = distribution_total + ...
            numel(distribution_groups) * numel(distribution_fields);
    end
end
[details, earned] = scoreCount(details, ...
    'residual_distribution_statistics', 2.0, distribution_passed, ...
    distribution_total, ...
    'Normal fit mean/std and direct residual RMSE/count were compared.');
score = score + earned;

report = finalizeReport(score, details, report_path);
end

function [passed, total] = compareMetricFields(actual, expected, group, ...
    fields, absolute_tolerances, relative_tolerance)
passed = 0;
total = numel(fields);
if ~isfield(actual, group) || ~isfield(expected, group)
    return;
end
actual_group = actual.(group);
expected_group = expected.(group);
for field_index = 1:numel(fields)
    field_name = fields{field_index};
    if numericFieldsEqual(actual_group, expected_group, field_name, ...
            absolute_tolerances(field_index), relative_tolerance)
        passed = passed + 1;
    end
end
end

function result = numericFieldsEqual(actual, expected, field_name, ...
    absolute_tolerance, relative_tolerance)
result = false;
if ~isfield(actual, field_name) || ~isfield(expected, field_name)
    return;
end
result = numericClose(actual.(field_name), expected.(field_name), ...
    absolute_tolerance, relative_tolerance);
end

function result = numericArrayFieldsEqual(actual, expected, field_name, ...
    absolute_tolerance, relative_tolerance)
result = false;
if ~isfield(actual, field_name) || ~isfield(expected, field_name)
    return;
end
actual_value = actual.(field_name);
expected_value = expected.(field_name);
if ~isnumeric(actual_value) || ~isnumeric(expected_value) || ...
        numel(actual_value) ~= numel(expected_value) || ...
        any(~isfinite(actual_value(:))) || any(~isfinite(expected_value(:)))
    return;
end
actual_value = double(actual_value(:));
expected_value = double(expected_value(:));
tolerance = max(absolute_tolerance, ...
    relative_tolerance * max(abs(expected_value), 1e-12));
result = all(abs(actual_value - expected_value) <= tolerance);
end

function result = numericClose(actual, expected, absolute_tolerance, ...
    relative_tolerance)
result = isnumeric(actual) && isnumeric(expected) && ...
    isscalar(actual) && isscalar(expected) && ...
    isfinite(actual) && isfinite(expected) && ...
    abs(double(actual) - double(expected)) <= max(absolute_tolerance, ...
    relative_tolerance * max(abs(double(expected)), 1e-12));
end

function result = booleanFieldsEqual(actual, expected, field_name)
result = false;
if ~isfield(actual, field_name) || ~isfield(expected, field_name)
    return;
end
actual_value = actual.(field_name);
expected_value = expected.(field_name);
result = (islogical(actual_value) || isnumeric(actual_value)) && ...
    (islogical(expected_value) || isnumeric(expected_value)) && ...
    isscalar(actual_value) && isscalar(expected_value) && ...
    logical(actual_value) == logical(expected_value);
end

function result = textFieldEquals(value, field_name, expected)
result = isfield(value, field_name) && ...
    isscalar(string(value.(field_name))) && ...
    string(value.(field_name)) == string(expected);
end

function result = trueField(value, field_name)
result = isfield(value, field_name) && ...
    (islogical(value.(field_name)) || isnumeric(value.(field_name))) && ...
    isscalar(value.(field_name)) && logical(value.(field_name));
end

function value = getFieldOrEmpty(container, field_name)
if isfield(container, field_name)
    value = container.(field_name);
else
    value = {};
end
end

function result = stringListsEqual(actual, expected)
try
    actual_strings = string(actual(:));
    expected_strings = string(expected(:));
    result = isequal(actual_strings, expected_strings);
catch
    result = false;
end
end

function [details, earned] = scoreChecks(details, criterion, ...
    max_points, checks, message)
[details, earned] = scoreCount(details, criterion, max_points, ...
    sum(checks), numel(checks), message);
end

function [details, earned] = scoreCount(details, criterion, max_points, ...
    passed, total, message)
if total <= 0
    earned = 0.0;
else
    earned = max_points * passed / total;
end
details = appendDetail(details, criterion, max_points, earned, ...
    passed, total, message);
end

function details = appendDetail(details, criterion, max_points, earned, ...
    passed, total, message)
item = struct;
item.criterion = criterion;
item.max_points = max_points;
item.earned_points = round(max(0, min(max_points, earned)), 6);
item.passed_checks = passed;
item.total_checks = total;
item.message = message;
details(end + 1, 1) = item;
end

function report = finalizeReport(score, details, report_path)
report = struct;
report.schema_version = '1.0';
report.record_id = 'matlab_nep_dft_energy_force_stress_validation';
report.score = round(max(0, min(100, score)), 6);
report.max_score = 100.0;
report.pass_score = 80.0;
report.passed = report.score + 1e-9 >= report.pass_score;
report.details = details;
report.evaluator = 'harness/evaluate_submission.m';

report_parent = fileparts(report_path);
if ~isempty(report_parent) && ~isfolder(report_parent)
    mkdir(report_parent);
end
writeJson(report_path, report);
end

function writeJson(path, value)
text = jsonencode(value, 'PrettyPrint', true);
fid = fopen(path, 'wt');
assert(fid >= 0, 'Cannot open JSON report: %s', path);
file_cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', text);
clear file_cleanup;
end
