function solution(input_dir, output_dir)
if nargin ~= 2
    error('solution requires exactly two inputs: input_dir and output_dir.');
end

input_dir = char(input_dir);
output_dir = char(output_dir);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

reference_file = fullfile(input_dir, 'test_reference.extxyz');
prediction_file = fullfile(input_dir, 'test_nep_predictions.h5');
training_file = fullfile(input_dir, 'training_history.csv');

reference = read_reference_extxyz(reference_file);
predictions = read_prediction_h5(prediction_file);
dataset = align_predictions(reference, predictions);
training = read_training_history(training_file);

prediction_before = baseline_prediction(dataset);
metrics_before = compute_metrics(dataset, prediction_before);
residual_rows = compute_residual_rows(dataset, prediction_before);
residual_stats_before = compute_residual_statistics(dataset, prediction_before);

metric_quality_pass = passes_metric_thresholds(metrics_before.overall);
training_assessment = build_training_assessment(training, metric_quality_pass);

calibration = solve_calibration(dataset);
prediction_after = apply_calibration(dataset, calibration.parameters);
metrics_after = compute_metrics(dataset, prediction_after);
residual_stats_after = compute_residual_statistics(dataset, prediction_after);

write_metrics_json(fullfile(output_dir, 'metrics.json'), metrics_before);
write_residuals_csv(fullfile(output_dir, 'residuals.csv'), residual_rows);
write_json_file(fullfile(output_dir, 'training_assessment.json'), training_assessment);
write_calibration_json( ...
    fullfile(output_dir, 'calibration.json'), ...
    calibration, ...
    metrics_before, ...
    metrics_after, ...
    residual_stats_before, ...
    residual_stats_after);
write_matlab_trace( ...
    fullfile(output_dir, 'matlab_trace.txt'), ...
    input_dir, ...
    output_dir, ...
    dataset, ...
    metrics_before, ...
    metrics_after, ...
    calibration, ...
    training_assessment);
end

function reference = read_reference_extxyz(file_path)
fid = fopen(file_path, 'r');
if fid < 0
    error('Failed to open reference file: %s', file_path);
end
cleanup_obj = onCleanup(@() fclose(fid));

allowed_materials = {'CsSnBr3', 'Cs2SnBr6', 'hetero_N2'};
structure_id = {};
material = {};
n_atoms = [];
energy_ref_eV = [];
stress_ref_GPa = [];
atom_id = {};
forces_ref_eV_A = {};

while true
    first_line = fgetl(fid);
    if ~ischar(first_line)
        break;
    end
    first_line = strtrim(first_line);
    if isempty(first_line)
        error('Unexpected blank line in reference file.');
    end

    atom_count = str2double(first_line);
    if ~isfinite(atom_count) || atom_count ~= round(atom_count) || atom_count <= 0
        error('Invalid atom count line in reference file: %s', first_line);
    end
    atom_count = double(atom_count);

    meta_line = fgetl(fid);
    if ~ischar(meta_line)
        error('Unexpected end of file while reading frame metadata.');
    end

    sid = extract_required_token(meta_line, '(?:^|\s)structure_id=([^\s]+)', 'structure_id');
    mat = extract_required_token(meta_line, '(?:^|\s)material=([^\s]+)', 'material');
    if ~ismember(mat, allowed_materials)
        error('Unsupported material label in reference file: %s', mat);
    end

    energy_token = extract_required_token(meta_line, '(?:^|\s)energy_ref_eV=([^\s]+)', 'energy_ref_eV');
    energy_value = str2double(energy_token);
    if ~isfinite(energy_value)
        error('Invalid energy_ref_eV value for structure %s.', sid);
    end

    stress_token = extract_optional_token(meta_line, 'stress_ref_GPa="([^"]+)"');
    if isempty(stress_token)
        stress_token = extract_required_token(meta_line, '(?:^|\s)stress_ref_GPa=([^\s]+)', 'stress_ref_GPa');
    end
    stress_value = sscanf(stress_token, '%f').';
    if numel(stress_value) ~= 6
        error('Expected 6 stress components for structure %s.', sid);
    end

    local_atom_id = zeros(atom_count, 1);
    local_forces = zeros(atom_count, 3);
    for atom_index = 1:atom_count
        atom_line = fgetl(fid);
        if ~ischar(atom_line)
            error('Unexpected end of file while reading atoms for structure %s.', sid);
        end
        tokens = strsplit(strtrim(atom_line));
        if numel(tokens) < 8
            error('Malformed atom line in structure %s: %s', sid, atom_line);
        end

        local_atom_id(atom_index, 1) = parse_scalar(tokens{5}, 'atom_id');
        local_forces(atom_index, 1) = parse_scalar(tokens{6}, 'fx_ref');
        local_forces(atom_index, 2) = parse_scalar(tokens{7}, 'fy_ref');
        local_forces(atom_index, 3) = parse_scalar(tokens{8}, 'fz_ref');
    end

    if numel(unique(local_atom_id)) ~= atom_count
        error('Duplicate atom_id detected in reference structure %s.', sid);
    end

    structure_id{end + 1, 1} = sid; %#ok<AGROW>
    material{end + 1, 1} = mat; %#ok<AGROW>
    n_atoms(end + 1, 1) = atom_count; %#ok<AGROW>
    energy_ref_eV(end + 1, 1) = energy_value; %#ok<AGROW>
    stress_ref_GPa(end + 1, :) = stress_value; %#ok<AGROW>
    atom_id{end + 1, 1} = local_atom_id; %#ok<AGROW>
    forces_ref_eV_A{end + 1, 1} = local_forces; %#ok<AGROW>
end

clear cleanup_obj;

if isempty(structure_id)
    error('Reference file contains no structures.');
end
if numel(unique(structure_id)) ~= numel(structure_id)
    error('Duplicate structure_id detected in reference file.');
end

reference = struct();
reference.structure_id = structure_id;
reference.material = material;
reference.n_atoms = n_atoms;
reference.energy_ref_eV = energy_ref_eV;
reference.energy_ref_eV_atom = energy_ref_eV ./ n_atoms;
reference.stress_ref_GPa = stress_ref_GPa;
reference.atom_id = atom_id;
reference.forces_ref_eV_A = forces_ref_eV_A;
end

function predictions = read_prediction_h5(file_path)
energy_pred_eV = read_h5_numeric_vector(file_path, '/energy_pred_eV', []);
n_structures = numel(energy_pred_eV);
structure_id = read_h5_text_vector(file_path, '/structure_id', n_structures);
stress_pred_GPa = read_h5_numeric_matrix(file_path, '/stress_pred_GPa', n_structures, 6);
offsets = read_h5_numeric_vector(file_path, '/offsets', n_structures + 1);
atom_id = read_h5_numeric_vector(file_path, '/atom_id', []);
n_atoms_total = numel(atom_id);
forces_pred_eV_A = read_h5_numeric_matrix(file_path, '/forces_pred_eV_A', n_atoms_total, 3);

if numel(unique(structure_id)) ~= numel(structure_id)
    error('Duplicate structure_id detected in prediction file.');
end
validate_integer_like(offsets, 'offsets');
validate_integer_like(atom_id, 'atom_id');
if offsets(1) ~= 0
    error('Prediction offsets must start at 0.');
end
if any(diff(offsets) < 0)
    error('Prediction offsets must be nondecreasing.');
end
if offsets(end) ~= n_atoms_total
    error('Prediction offsets(end) must equal numel(atom_id).');
end

predictions = struct();
predictions.structure_id = structure_id;
predictions.energy_pred_eV = energy_pred_eV;
predictions.stress_pred_GPa = stress_pred_GPa;
predictions.offsets = offsets;
predictions.atom_id = atom_id;
predictions.forces_pred_eV_A = forces_pred_eV_A;
end

function dataset = align_predictions(reference, predictions)
n_structures = numel(reference.structure_id);
if numel(predictions.structure_id) ~= n_structures
    error('Reference and prediction files contain different numbers of structures.');
end

if ~isequal(sort(reference.structure_id), sort(predictions.structure_id))
    error('Reference and prediction structure_id sets do not match exactly.');
end

prediction_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for idx = 1:n_structures
    prediction_map(predictions.structure_id{idx}) = idx;
end

aligned_energy = zeros(n_structures, 1);
aligned_stress = zeros(n_structures, 6);
aligned_forces = cell(n_structures, 1);

for idx = 1:n_structures
    sid = reference.structure_id{idx};
    pred_idx = prediction_map(sid);

    start_idx = predictions.offsets(pred_idx) + 1;
    end_idx = predictions.offsets(pred_idx + 1);
    atom_count = end_idx - start_idx + 1;
    if atom_count ~= reference.n_atoms(idx)
        error('Atom count mismatch for structure %s.', sid);
    end

    pred_atom_id = predictions.atom_id(start_idx:end_idx);
    pred_forces = predictions.forces_pred_eV_A(start_idx:end_idx, :);
    ref_atom_id = reference.atom_id{idx};

    if numel(unique(pred_atom_id)) ~= numel(pred_atom_id)
        error('Duplicate predicted atom_id detected for structure %s.', sid);
    end
    if numel(pred_atom_id) ~= numel(ref_atom_id) || ...
            ~isequal(sort(pred_atom_id(:)), sort(ref_atom_id(:)))
        error('Predicted atom_id set mismatch for structure %s.', sid);
    end

    [is_present, locations] = ismember(ref_atom_id(:), pred_atom_id(:));
    if ~all(is_present) || any(locations == 0)
        error('Failed to align predicted forces by atom_id for structure %s.', sid);
    end

    aligned_energy(idx, 1) = predictions.energy_pred_eV(pred_idx);
    aligned_stress(idx, :) = predictions.stress_pred_GPa(pred_idx, :);
    aligned_forces{idx, 1} = pred_forces(locations, :);
end

dataset = reference;
dataset.energy_pred_eV = aligned_energy;
dataset.energy_pred_eV_atom = aligned_energy ./ dataset.n_atoms;
dataset.stress_pred_GPa = aligned_stress;
dataset.forces_pred_eV_A = aligned_forces;
end

function training = read_training_history(file_path)
training_table = readtable(file_path);
required_columns = { ...
    'step', ...
    'train_total_loss', ...
    'test_total_loss', ...
    'train_energy_loss', ...
    'test_energy_loss', ...
    'train_force_loss', ...
    'test_force_loss', ...
    'train_stress_loss', ...
    'test_stress_loss'};

missing_columns = required_columns(~ismember(required_columns, training_table.Properties.VariableNames));
if ~isempty(missing_columns)
    error('Missing required training history columns: %s', strjoin(missing_columns, ', '));
end

training = struct();
training.n_records = height(training_table);
training.step = double(training_table.step(:));
training.train_total_loss = double(training_table.train_total_loss(:));
training.test_total_loss = double(training_table.test_total_loss(:));
training.train_energy_loss = double(training_table.train_energy_loss(:));
training.test_energy_loss = double(training_table.test_energy_loss(:));
training.train_force_loss = double(training_table.train_force_loss(:));
training.test_force_loss = double(training_table.test_force_loss(:));
training.train_stress_loss = double(training_table.train_stress_loss(:));
training.test_stress_loss = double(training_table.test_stress_loss(:));
end

function prediction = baseline_prediction(dataset)
prediction = struct();
prediction.energy_pred_eV = dataset.energy_pred_eV;
prediction.energy_pred_eV_atom = dataset.energy_pred_eV_atom;
prediction.stress_pred_GPa = dataset.stress_pred_GPa;
prediction.forces_pred_eV_A = dataset.forces_pred_eV_A;
end

function metrics = compute_metrics(dataset, prediction)
group_names = {'overall', 'CsSnBr3', 'Cs2SnBr6', 'hetero_N2'};
metrics = struct();

for idx = 1:numel(group_names)
    group_name = group_names{idx};
    if strcmp(group_name, 'overall')
        mask = true(numel(dataset.structure_id), 1);
    else
        mask = strcmp(dataset.material, group_name);
    end
    group_indices = find(mask);
    metrics.(group_name) = compute_group_metrics(dataset, prediction, group_indices);
end
end

function group_metrics = compute_group_metrics(dataset, prediction, indices)
if isempty(indices)
    error('Encountered an empty metric group.');
end

ref_energy = dataset.energy_ref_eV_atom(indices);
pred_energy = prediction.energy_pred_eV_atom(indices);
energy_residual_meV_atom = (pred_energy - ref_energy) * 1000.0;

ref_force = flatten_force_cells(dataset.forces_ref_eV_A(indices));
pred_force = flatten_force_cells(prediction.forces_pred_eV_A(indices));

ref_stress = dataset.stress_ref_GPa(indices, :);
pred_stress = prediction.stress_pred_GPa(indices, :);

group_metrics = struct();
group_metrics.count = numel(indices);
group_metrics.energy_RMSE_meV_atom = sqrt(mean(energy_residual_meV_atom .^ 2));
group_metrics.energy_R2 = compute_r2(ref_energy, pred_energy);
group_metrics.force_RMSE_eV_A = sqrt(mean((pred_force - ref_force) .^ 2));
group_metrics.force_R2 = compute_r2(ref_force, pred_force);
group_metrics.stress_RMSE_GPa = sqrt(mean((pred_stress(:) - ref_stress(:)) .^ 2));
group_metrics.stress_R2 = compute_r2(ref_stress(:), pred_stress(:));
end

function residual_rows = compute_residual_rows(dataset, prediction)
n_structures = numel(dataset.structure_id);
residual_rows = repmat(struct( ...
    'structure_id', '', ...
    'material', '', ...
    'n_atoms', 0, ...
    'energy_ref_eV', 0, ...
    'energy_pred_eV', 0, ...
    'energy_error_meV_atom', 0, ...
    'force_rmse_eV_A', 0, ...
    'stress_rmse_GPa', 0, ...
    'max_force_abs_error_eV_A', 0), n_structures, 1);

for idx = 1:n_structures
    force_diff = prediction.forces_pred_eV_A{idx} - dataset.forces_ref_eV_A{idx};
    stress_diff = prediction.stress_pred_GPa(idx, :) - dataset.stress_ref_GPa(idx, :);

    residual_rows(idx).structure_id = dataset.structure_id{idx};
    residual_rows(idx).material = dataset.material{idx};
    residual_rows(idx).n_atoms = dataset.n_atoms(idx);
    residual_rows(idx).energy_ref_eV = dataset.energy_ref_eV(idx);
    residual_rows(idx).energy_pred_eV = prediction.energy_pred_eV(idx);
    residual_rows(idx).energy_error_meV_atom = ...
        (prediction.energy_pred_eV(idx) - dataset.energy_ref_eV(idx)) / dataset.n_atoms(idx) * 1000.0;
    residual_rows(idx).force_rmse_eV_A = sqrt(mean(force_diff(:) .^ 2));
    residual_rows(idx).stress_rmse_GPa = sqrt(mean(stress_diff(:) .^ 2));
    residual_rows(idx).max_force_abs_error_eV_A = max(abs(force_diff(:)));
end
end

function residual_stats = compute_residual_statistics(dataset, prediction)
[energy_residual, force_residual, stress_residual] = collect_residual_vectors(dataset, prediction);

residual_stats = struct();
residual_stats.energy_meV_atom = fit_normal_statistics(energy_residual);
residual_stats.force_eV_A = fit_normal_statistics(force_residual);
residual_stats.stress_GPa = fit_normal_statistics(stress_residual);
end

function [energy_residual, force_residual, stress_residual] = collect_residual_vectors(dataset, prediction)
energy_residual = (prediction.energy_pred_eV - dataset.energy_ref_eV) ./ dataset.n_atoms * 1000.0;

force_diff = cell(numel(dataset.forces_ref_eV_A), 1);
for idx = 1:numel(force_diff)
    force_diff{idx} = prediction.forces_pred_eV_A{idx} - dataset.forces_ref_eV_A{idx};
end
force_residual = flatten_force_cells(force_diff);
stress_residual = prediction.stress_pred_GPa(:) - dataset.stress_ref_GPa(:);
end

function stats = fit_normal_statistics(residual)
residual = residual(:);
distribution = fitdist(residual, 'Normal');

stats = struct();
stats.mean = distribution.mu;
stats.std = distribution.sigma;
stats.rmse = sqrt(mean(residual .^ 2));
stats.count = numel(residual);
end

function pass = passes_metric_thresholds(overall_metrics)
pass = ...
    overall_metrics.energy_RMSE_meV_atom <= 2.0 && ...
    overall_metrics.force_RMSE_eV_A <= 0.08 && ...
    overall_metrics.stress_RMSE_GPa <= 0.10 && ...
    overall_metrics.energy_R2 >= 0.99 && ...
    overall_metrics.force_R2 >= 0.70 && ...
    overall_metrics.stress_R2 >= 0.95;
end

function assessment = build_training_assessment(training, metric_quality_pass)
initial_train = training.train_total_loss(1);
final_train = training.train_total_loss(end);
initial_test = training.test_total_loss(1);
final_test = training.test_total_loss(end);
[best_test, best_idx] = min(training.test_total_loss);

train_reduction = compute_reduction_fraction(initial_train, final_train);
test_reduction = compute_reduction_fraction(initial_test, final_test);
final_to_best = safe_ratio(final_test, best_test);
final_test_to_train = safe_ratio(final_test, final_train);

converged = train_reduction >= 0.80 && test_reduction >= 0.80;
no_overfitting = final_to_best <= 1.10 && final_test_to_train <= 1.35;

assessment = struct();
assessment.schema_version = '1.0';
assessment.synthetic_data = true;
assessment.n_records = training.n_records;
assessment.first_step = training.step(1);
assessment.last_step = training.step(end);
assessment.initial_train_total_loss = initial_train;
assessment.final_train_total_loss = final_train;
assessment.initial_test_total_loss = initial_test;
assessment.final_test_total_loss = final_test;
assessment.train_loss_reduction_fraction = train_reduction;
assessment.test_loss_reduction_fraction = test_reduction;
assessment.best_test_step = training.step(best_idx);
assessment.final_to_best_test_loss_ratio = final_to_best;
assessment.final_test_to_train_loss_ratio = final_test_to_train;
assessment.converged = converged;
assessment.no_overfitting = no_overfitting;
assessment.metric_quality_pass = metric_quality_pass;
assessment.reliable_for_thermal_transport = converged && no_overfitting && metric_quality_pass;
end

function calibration = solve_calibration(dataset)
energy_ref = dataset.energy_ref_eV_atom(:);
energy_pred = dataset.energy_pred_eV_atom(:);
force_ref = flatten_force_cells(dataset.forces_ref_eV_A);
force_pred = flatten_force_cells(dataset.forces_pred_eV_A);
stress_ref = dataset.stress_ref_GPa(:);
stress_pred = dataset.stress_pred_GPa(:);

n_energy = numel(energy_ref);
n_force = numel(force_ref);
n_stress = numel(stress_ref);

energy_scale = 0.01;
force_scale = 0.10;
stress_scale = 0.10;
lambda = 0.01;

theta0 = [0.0; 1.0; 1.0];
theta_center = [0.0; 1.0; 1.0];
theta_scales = [0.005; 0.05; 0.05];
lb = [-0.005; 0.80; 0.80];
ub = [0.005; 1.20; 1.20];

energy_weight = 1.0 / (energy_scale * sqrt(n_energy));
force_weight = 1.0 / (force_scale * sqrt(n_force));
stress_weight = 1.0 / (stress_scale * sqrt(n_stress));
sqrt_lambda = sqrt(lambda);

design_energy = [energy_weight * ones(n_energy, 1), zeros(n_energy, 2)];
target_energy = energy_weight * (energy_ref - energy_pred);

design_force = [zeros(n_force, 1), force_weight * force_pred, zeros(n_force, 1)];
target_force = force_weight * force_ref;

design_stress = [zeros(n_stress, 2), stress_weight * stress_pred];
target_stress = stress_weight * stress_ref;

design_regularization = diag(sqrt_lambda ./ theta_scales);
target_regularization = design_regularization * theta_center;

design_matrix = [design_energy; design_force; design_stress; design_regularization];
target_vector = [target_energy; target_force; target_stress; target_regularization];

options = optimoptions('lsqlin', 'Display', 'off');
strictly_convex = rank(design_matrix) == size(design_matrix, 2);
solver_name = 'lsqlin';
theta = theta0;
resnorm = NaN;
exitflag = NaN;
output = struct();

try
    [theta, resnorm, ~, exitflag, output] = lsqlin( ...
        design_matrix, ...
        target_vector, ...
        [], [], [], [], ...
        lb, ...
        ub, ...
        theta0, ...
        options);
catch solver_error
    error('Calibration solver failed: %s', solver_error.message);
end

iterations = get_struct_field_or_default(output, 'iterations', 0);

calibration = struct();
calibration.schema_version = '1.0';
calibration.synthetic_data = true;
calibration.method = 'bounded_regularized_lsqlin';
calibration.objective = struct( ...
    'energy_denominator_eV_atom', energy_scale, ...
    'force_denominator_eV_A', force_scale, ...
    'stress_denominator_GPa', stress_scale, ...
    'equal_block_weighting', true, ...
    'regularization_lambda', lambda, ...
    'regularization_center', struct('b_E', theta_center(1), 's_F', theta_center(2), 's_S', theta_center(3)), ...
    'regularization_scales', struct('b_E', theta_scales(1), 's_F', theta_scales(2), 's_S', theta_scales(3)));
calibration.bounds = struct( ...
    'b_E', struct('lower', lb(1), 'upper', ub(1)), ...
    's_F', struct('lower', lb(2), 'upper', ub(2)), ...
    's_S', struct('lower', lb(3), 'upper', ub(3)));
calibration.parameters = struct( ...
    'energy_bias_eV_atom', theta(1), ...
    'force_scale', theta(2), ...
    'stress_scale', theta(3));
calibration.solver = struct( ...
    'name', solver_name, ...
    'exitflag', exitflag, ...
    'iterations', iterations, ...
    'resnorm', resnorm);
calibration.strictly_convex = strictly_convex;
calibration.design_matrix_rank = rank(design_matrix);
calibration.objective_value = resnorm;
calibration.unique_solution = strictly_convex && exitflag > 0;
end

function prediction = apply_calibration(dataset, parameters)
prediction = struct();
prediction.energy_pred_eV_atom = dataset.energy_pred_eV_atom + parameters.energy_bias_eV_atom;
prediction.energy_pred_eV = prediction.energy_pred_eV_atom .* dataset.n_atoms;
prediction.stress_pred_GPa = parameters.stress_scale * dataset.stress_pred_GPa;
prediction.forces_pred_eV_A = cellfun( ...
    @(force_matrix) parameters.force_scale * force_matrix, ...
    dataset.forces_pred_eV_A, ...
    'UniformOutput', false);
end

function write_metrics_json(file_path, metrics)
payload = struct();
payload.schema_version = '1.0';
payload.synthetic_data = true;
payload.stress_voigt_order = {'xx', 'yy', 'zz', 'yz', 'xz', 'xy'};
payload.overall = metrics.overall;
payload.CsSnBr3 = metrics.CsSnBr3;
payload.Cs2SnBr6 = metrics.Cs2SnBr6;
payload.hetero_N2 = metrics.hetero_N2;
write_json_file(file_path, payload);
end

function write_calibration_json(file_path, calibration, metrics_before, metrics_after, residual_stats_before, residual_stats_after)
payload = calibration;
payload.metrics_before = struct( ...
    'overall', metrics_before.overall, ...
    'CsSnBr3', metrics_before.CsSnBr3, ...
    'Cs2SnBr6', metrics_before.Cs2SnBr6, ...
    'hetero_N2', metrics_before.hetero_N2);
payload.metrics_after = struct( ...
    'overall', metrics_after.overall, ...
    'CsSnBr3', metrics_after.CsSnBr3, ...
    'Cs2SnBr6', metrics_after.Cs2SnBr6, ...
    'hetero_N2', metrics_after.hetero_N2);
payload.residual_statistics_before = residual_stats_before;
payload.residual_statistics_after = residual_stats_after;
write_json_file(file_path, payload);
end

function write_json_file(file_path, payload)
fid = fopen(file_path, 'w');
if fid < 0
    error('Failed to open JSON output file for writing: %s', file_path);
end
cleanup_obj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(payload));
clear cleanup_obj;
end

function write_residuals_csv(file_path, residual_rows)
fid = fopen(file_path, 'w');
if fid < 0
    error('Failed to open residual CSV for writing: %s', file_path);
end
cleanup_obj = onCleanup(@() fclose(fid));

fprintf(fid, 'structure_id,material,n_atoms,energy_ref_eV,energy_pred_eV,energy_error_meV_atom,force_rmse_eV_A,stress_rmse_GPa,max_force_abs_error_eV_A\n');
for idx = 1:numel(residual_rows)
    row = residual_rows(idx);
    fprintf(fid, '%s,%s,%d,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n', ...
        row.structure_id, ...
        row.material, ...
        row.n_atoms, ...
        row.energy_ref_eV, ...
        row.energy_pred_eV, ...
        row.energy_error_meV_atom, ...
        row.force_rmse_eV_A, ...
        row.stress_rmse_GPa, ...
        row.max_force_abs_error_eV_A);
end
clear cleanup_obj;
end

function write_matlab_trace(file_path, input_dir, output_dir, dataset, metrics_before, metrics_after, calibration, training_assessment)
toolbox_versions = get_toolbox_versions();
lines = { ...
    'synthetic_data=true', ...
    sprintf('input_dir=%s', input_dir), ...
    sprintf('output_dir=%s', output_dir), ...
    sprintf('matlab_version=%s', version), ...
    sprintf('matlab_release=%s', version('-release')), ...
    sprintf('optimization_toolbox_version=%s', toolbox_versions.optimization), ...
    sprintf('statistics_toolbox_version=%s', toolbox_versions.statistics), ...
    'toolbox_calls=h5read,readtable,lsqlin,fitdist,jsonencode', ...
    sprintf('n_structures=%d', numel(dataset.structure_id)), ...
    sprintf('total_atoms=%d', sum(dataset.n_atoms)), ...
    sprintf('overall_energy_rmse_before_meV_atom=%.15g', metrics_before.overall.energy_RMSE_meV_atom), ...
    sprintf('overall_force_rmse_before_eV_A=%.15g', metrics_before.overall.force_RMSE_eV_A), ...
    sprintf('overall_stress_rmse_before_GPa=%.15g', metrics_before.overall.stress_RMSE_GPa), ...
    sprintf('overall_energy_rmse_after_meV_atom=%.15g', metrics_after.overall.energy_RMSE_meV_atom), ...
    sprintf('overall_force_rmse_after_eV_A=%.15g', metrics_after.overall.force_RMSE_eV_A), ...
    sprintf('overall_stress_rmse_after_GPa=%.15g', metrics_after.overall.stress_RMSE_GPa), ...
    sprintf('calibration_energy_bias_eV_atom=%.15g', calibration.parameters.energy_bias_eV_atom), ...
    sprintf('calibration_force_scale=%.15g', calibration.parameters.force_scale), ...
    sprintf('calibration_stress_scale=%.15g', calibration.parameters.stress_scale), ...
    sprintf('calibration_exitflag=%.15g', calibration.solver.exitflag), ...
    sprintf('calibration_iterations=%.15g', calibration.solver.iterations), ...
    sprintf('calibration_resnorm=%.15g', calibration.solver.resnorm), ...
    sprintf('calibration_unique_solution=%d', calibration.unique_solution), ...
    sprintf('metric_quality_pass=%d', training_assessment.metric_quality_pass), ...
    sprintf('reliable_for_thermal_transport=%d', training_assessment.reliable_for_thermal_transport), ...
    'files_written=metrics.json,residuals.csv,training_assessment.json,calibration.json,matlab_trace.txt'};

fid = fopen(file_path, 'w');
if fid < 0
    error('Failed to open matlab_trace.txt for writing: %s', file_path);
end
cleanup_obj = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', lines{:});
clear cleanup_obj;
end

function toolbox_versions = get_toolbox_versions()
installed = ver;
toolbox_versions = struct();
toolbox_versions.optimization = lookup_toolbox_version(installed, 'Optimization Toolbox');
toolbox_versions.statistics = lookup_toolbox_version(installed, 'Statistics and Machine Learning Toolbox');
end

function version_string = lookup_toolbox_version(installed, toolbox_name)
names = {installed.Name};
match_idx = find(strcmp(names, toolbox_name), 1);
if isempty(match_idx)
    version_string = 'not_found';
else
    version_string = installed(match_idx).Version;
end
end

function vector = read_h5_numeric_vector(file_path, dataset_name, expected_length)
raw = h5read(file_path, dataset_name);
vector = squeeze(double(raw));
vector = vector(:);
if ~isempty(expected_length) && numel(vector) ~= expected_length
    error('Dataset %s has length %d; expected %d.', dataset_name, numel(vector), expected_length);
end
end

function matrix = read_h5_numeric_matrix(file_path, dataset_name, expected_rows, expected_cols)
raw = h5read(file_path, dataset_name);
matrix = squeeze(double(raw));
matrix_size = size(matrix);

if isvector(matrix)
    if expected_rows == 1 || expected_cols == 1
        matrix = reshape(matrix, expected_rows, expected_cols);
    else
        error('Dataset %s is unexpectedly one-dimensional.', dataset_name);
    end
elseif isequal(matrix_size, [expected_rows, expected_cols])
    % Already normalized.
elseif isequal(matrix_size, [expected_cols, expected_rows])
    matrix = matrix.';
else
    error('Dataset %s has shape [%s]; expected [%d %d] or [%d %d].', ...
        dataset_name, num2str(matrix_size), expected_rows, expected_cols, expected_cols, expected_rows);
end
end

function text = read_h5_text_vector(file_path, dataset_name, expected_count)
raw = h5read(file_path, dataset_name);
text = normalize_text_vector(raw, expected_count, dataset_name);
end

function text = normalize_text_vector(raw, expected_count, dataset_name)
if isstring(raw)
    text = cellstr(raw(:));
elseif iscell(raw)
    text = cell(numel(raw), 1);
    for idx = 1:numel(raw)
        item = raw{idx};
        if isstring(item)
            item = char(item);
        elseif isnumeric(item) || islogical(item)
            item = char(item(:).');
        elseif ~ischar(item)
            error('Unsupported HDF5 text cell element type in %s.', dataset_name);
        end
        text{idx, 1} = clean_text(item);
    end
elseif ischar(raw)
    text = char_matrix_to_cell(raw, expected_count, dataset_name);
elseif isnumeric(raw) || islogical(raw)
    if isvector(raw)
        if expected_count ~= 1
            error('Cannot interpret numeric vector dataset %s as %d strings.', dataset_name, expected_count);
        end
        text = {clean_text(char(raw(:).'))};
    else
        text = char_matrix_to_cell(char(raw), expected_count, dataset_name);
    end
else
    error('Unsupported HDF5 text dataset type for %s.', dataset_name);
end

text = text(:);
for idx = 1:numel(text)
    text{idx} = clean_text(text{idx});
end

if numel(text) ~= expected_count
    error('Dataset %s produced %d strings; expected %d.', dataset_name, numel(text), expected_count);
end
end

function text = char_matrix_to_cell(raw, expected_count, dataset_name)
raw(raw == 0) = ' ';
if size(raw, 1) == expected_count
    text = cellstr(raw);
elseif size(raw, 2) == expected_count
    text = cellstr(raw.');
elseif expected_count == 1
    text = {clean_text(raw(:).')};
else
    error('Cannot normalize character array dataset %s to %d strings.', dataset_name, expected_count);
end
text = text(:);
end

function cleaned = clean_text(value)
cleaned = char(value);
cleaned = cleaned(:).';
cleaned(cleaned == 0) = ' ';
cleaned = strtrim(cleaned);
end

function values = flatten_force_cells(force_cells)
if isempty(force_cells)
    values = zeros(0, 1);
    return;
end
values = vertcat(force_cells{:});
values = values(:);
end

function value = compute_r2(reference, prediction)
reference = reference(:);
prediction = prediction(:);
residual = prediction - reference;
centered_reference = reference - mean(reference);
sst = sum(centered_reference .^ 2);
if sst == 0
    if all(abs(residual) <= eps(max(1.0, max(abs(reference)))))
        value = 1.0;
    else
        value = NaN;
    end
    return;
end
value = 1.0 - sum(residual .^ 2) / sst;
end

function fraction = compute_reduction_fraction(initial_value, final_value)
if initial_value == 0
    if final_value == 0
        fraction = 1.0;
    else
        fraction = -Inf;
    end
else
    fraction = (initial_value - final_value) / initial_value;
end
end

function ratio = safe_ratio(numerator, denominator)
if denominator == 0
    if numerator == 0
        ratio = 1.0;
    else
        ratio = Inf;
    end
else
    ratio = numerator / denominator;
end
end

function value = parse_scalar(token, field_name)
value = str2double(token);
if ~isfinite(value)
    error('Invalid numeric token for %s: %s', field_name, token);
end
value = double(value);
end

function token = extract_required_token(text, pattern, field_name)
token = extract_optional_token(text, pattern);
if isempty(token)
    error('Missing required %s in metadata line: %s', field_name, text);
end
end

function token = extract_optional_token(text, pattern)
match = regexp(text, pattern, 'tokens', 'once');
if isempty(match)
    token = '';
else
    token = match{1};
end
end

function validate_integer_like(values, field_name)
if any(abs(values - round(values)) > 1e-9)
    error('Expected integer-like values in %s.', field_name);
end
end

function value = get_struct_field_or_default(struct_value, field_name, default_value)
if isstruct(struct_value) && isfield(struct_value, field_name)
    value = struct_value.(field_name);
else
    value = default_value;
end
end
