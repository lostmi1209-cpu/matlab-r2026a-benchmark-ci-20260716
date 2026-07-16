function solution(input_dir, submission_dir)
%SOLUTION Evaluate aligned synthetic NEP predictions and calibrate them.

if nargin ~= 2
    error('solution:InvalidArguments', ...
        'Expected solution(input_dir, submission_dir).');
end
input_dir = char(input_dir);
submission_dir = char(submission_dir);
if ~isfolder(input_dir)
    error('solution:MissingInputDirectory', 'Input directory not found: %s', input_dir);
end
if ~isfolder(submission_dir)
    mkdir(submission_dir);
end

reference_path = fullfile(input_dir, 'test_reference.extxyz');
prediction_path = fullfile(input_dir, 'test_nep_predictions.h5');
history_path = fullfile(input_dir, 'training_history.csv');

fprintf('MATLAB version: %s\n', version);
fprintf('Reading and aligning public synthetic inputs.\n');
reference = read_reference_extxyz(reference_path);
prediction = read_prediction_h5(prediction_path);
data = align_predictions(reference, prediction);

metrics_before = compute_all_metrics(data, data.energy_pred_atom, ...
    data.force_pred, data.stress_pred);

metrics_document = struct();
metrics_document.schema_version = '1.0';
metrics_document.synthetic_data = true;
metrics_document.stress_voigt_order = {'xx', 'yy', 'zz', 'yz', 'xz', 'xy'};
metrics_document = append_group_metrics(metrics_document, metrics_before);
write_json(fullfile(submission_dir, 'metrics.json'), metrics_document);

write_residual_table(fullfile(submission_dir, 'residuals.csv'), data);

training = assess_training(history_path, metrics_before.overall);
write_json(fullfile(submission_dir, 'training_assessment.json'), training);

fprintf('Calling Optimization Toolbox lsqlin.\n');
[theta, solver_info, convexity] = calibrate_predictions(data);
energy_pred_after = data.energy_pred_atom + theta(1);
force_pred_after = scale_cell_matrices(data.force_pred, theta(2));
stress_pred_after = theta(3) .* data.stress_pred;
metrics_after = compute_all_metrics(data, energy_pred_after, ...
    force_pred_after, stress_pred_after);

fprintf('Calling Statistics and Machine Learning Toolbox fitdist.\n');
stats_before = residual_statistics(data, data.energy_pred_atom, ...
    data.force_pred, data.stress_pred);
stats_after = residual_statistics(data, energy_pred_after, ...
    force_pred_after, stress_pred_after);

calibration = build_calibration_document(theta, solver_info, convexity, ...
    metrics_before, metrics_after, stats_before, stats_after);
write_json(fullfile(submission_dir, 'calibration.json'), calibration);

fprintf(['Overall metrics before calibration: energy RMSE %.12g meV/atom, ' ...
    'force RMSE %.12g eV/A, stress RMSE %.12g GPa.\n'], ...
    metrics_before.overall.energy_RMSE_meV_atom, ...
    metrics_before.overall.force_RMSE_eV_A, ...
    metrics_before.overall.stress_RMSE_GPa);
fprintf('Calibration parameters: b_E=%.17g, s_F=%.17g, s_S=%.17g.\n', ...
    theta(1), theta(2), theta(3));
fprintf('lsqlin exitflag=%d, resnorm=%.17g.\n', ...
    solver_info.exitflag, solver_info.resnorm);
fprintf('Outputs written to: %s\n', submission_dir);
end

function frames = read_reference_extxyz(path)
fid = fopen(path, 'r');
if fid < 0
    error('solution:OpenFailed', 'Unable to open reference file: %s', path);
end
cleanup = onCleanup(@() fclose(fid));
frames = struct('structure_id', {}, 'material', {}, 'n_atoms', {}, ...
    'energy_ref_eV', {}, 'stress_ref_GPa', {}, 'atom_id', {}, ...
    'force_ref_eV_A', {});

while true
    count_line = fgetl(fid);
    while ischar(count_line) && isempty(strtrim(count_line))
        count_line = fgetl(fid);
    end
    if ~ischar(count_line)
        break;
    end
    n_atoms = str2double(strtrim(count_line));
    if ~isfinite(n_atoms) || n_atoms <= 0 || n_atoms ~= round(n_atoms)
        error('solution:InvalidExtxyz', 'Invalid atom count line: %s', count_line);
    end
    metadata = fgetl(fid);
    if ~ischar(metadata)
        error('solution:InvalidExtxyz', 'Unexpected EOF in frame metadata.');
    end

    structure_id = metadata_value(metadata, 'structure_id');
    material = metadata_value(metadata, 'material');
    energy_ref = str2double(metadata_value(metadata, 'energy_ref_eV'));
    stress_ref = sscanf(metadata_value(metadata, 'stress_ref_GPa'), '%f').';
    if ~isfinite(energy_ref) || numel(stress_ref) ~= 6 || any(~isfinite(stress_ref))
        error('solution:InvalidExtxyz', ...
            'Invalid energy or stress metadata for structure %s.', structure_id);
    end

    atom_id = zeros(n_atoms, 1);
    force_ref = zeros(n_atoms, 3);
    for atom_index = 1:n_atoms
        atom_line = fgetl(fid);
        if ~ischar(atom_line)
            error('solution:InvalidExtxyz', ...
                'Unexpected EOF in atoms for structure %s.', structure_id);
        end
        fields = strsplit(strtrim(atom_line));
        if numel(fields) ~= 8
            error('solution:InvalidExtxyz', ...
                'Expected 8 atom fields in structure %s.', structure_id);
        end
        numeric_fields = str2double(fields(2:8));
        if any(~isfinite(numeric_fields))
            error('solution:InvalidExtxyz', ...
                'Non-numeric atom data in structure %s.', structure_id);
        end
        current_id = numeric_fields(4);
        if current_id ~= round(current_id)
            error('solution:InvalidExtxyz', ...
                'Non-integer atom_id in structure %s.', structure_id);
        end
        atom_id(atom_index) = current_id;
        force_ref(atom_index, :) = numeric_fields(5:7);
    end
    if numel(unique(atom_id)) ~= n_atoms
        error('solution:DuplicateAtomId', ...
            'Duplicate reference atom_id in structure %s.', structure_id);
    end

    frame = struct();
    frame.structure_id = structure_id;
    frame.material = material;
    frame.n_atoms = n_atoms;
    frame.energy_ref_eV = energy_ref;
    frame.stress_ref_GPa = stress_ref;
    frame.atom_id = atom_id;
    frame.force_ref_eV_A = force_ref;
    frames(end + 1, 1) = frame; %#ok<AGROW>
end
clear cleanup;

if isempty(frames)
    error('solution:InvalidExtxyz', 'Reference file contains no structures.');
end
ids = string({frames.structure_id}).';
if numel(unique(ids)) ~= numel(ids)
    error('solution:DuplicateStructureId', ...
        'Reference structure_id values are not unique.');
end
end

function value = metadata_value(line, key)
escaped_key = regexptranslate('escape', key);
token = regexp(line, [escaped_key '="([^"]*)"'], 'tokens', 'once');
if ~isempty(token)
    value = token{1};
    return;
end
token = regexp(line, [escaped_key '=([^\s]+)'], 'tokens', 'once');
if isempty(token)
    error('solution:MissingMetadata', 'Missing metadata field %s.', key);
end
value = token{1};
end

function prediction = read_prediction_h5(path)
energy = double(h5read(path, '/energy_pred_eV'));
energy = energy(:);
n_structures = numel(energy);
if n_structures == 0 || any(~isfinite(energy))
    error('solution:InvalidH5', 'Invalid /energy_pred_eV dataset.');
end

structure_id = normalize_text_vector(h5read(path, '/structure_id'), ...
    n_structures, '/structure_id');
stress = normalize_component_matrix(h5read(path, '/stress_pred_GPa'), ...
    n_structures, 6, '/stress_pred_GPa');
offsets = double(h5read(path, '/offsets'));
offsets = offsets(:);
atom_id = double(h5read(path, '/atom_id'));
atom_id = atom_id(:);
n_atoms_total = numel(atom_id);
forces = normalize_component_matrix(h5read(path, '/forces_pred_eV_A'), ...
    n_atoms_total, 3, '/forces_pred_eV_A');

if numel(offsets) ~= n_structures + 1 || any(~isfinite(offsets)) || ...
        any(abs(offsets - round(offsets)) > 1e-9) || offsets(1) ~= 0 || ...
        any(diff(offsets) < 0) || offsets(end) ~= n_atoms_total
    error('solution:InvalidH5', 'Invalid zero-based /offsets dataset.');
end
if any(~isfinite(atom_id)) || any(abs(atom_id - round(atom_id)) > 1e-9)
    error('solution:InvalidH5', 'Invalid /atom_id dataset.');
end
if any(~isfinite(stress(:))) || any(~isfinite(forces(:)))
    error('solution:InvalidH5', 'Prediction arrays contain non-finite values.');
end
if numel(unique(structure_id)) ~= n_structures
    error('solution:DuplicateStructureId', ...
        'Prediction structure_id values are not unique.');
end

prediction = struct();
prediction.structure_id = structure_id;
prediction.energy_pred_eV = energy;
prediction.stress_pred_GPa = stress;
prediction.offsets = round(offsets);
prediction.atom_id = round(atom_id);
prediction.forces_pred_eV_A = forces;
end

function matrix = normalize_component_matrix(raw, n_rows, n_components, name)
matrix = squeeze(double(raw));
if isequal(size(matrix), [n_rows, n_components])
    return;
end
if isequal(size(matrix), [n_components, n_rows])
    matrix = matrix.';
    return;
end
if n_rows == 1 && numel(matrix) == n_components
    matrix = reshape(matrix, 1, n_components);
    return;
end
error('solution:InvalidH5Shape', ...
    '%s has shape incompatible with logical %d-by-%d.', ...
    name, n_rows, n_components);
end

function ids = normalize_text_vector(raw, expected_count, name)
if isstring(raw)
    ids = raw(:);
elseif iscell(raw)
    ids = strings(numel(raw), 1);
    for index = 1:numel(raw)
        ids(index) = decode_text_item(raw{index}, name);
    end
elseif ischar(raw) || isnumeric(raw)
    ids = decode_text_array(raw, expected_count, name);
else
    error('solution:InvalidH5String', 'Unsupported string type for %s.', name);
end
ids = strip(ids(:));
if numel(ids) ~= expected_count || any(strlength(ids) == 0)
    error('solution:InvalidH5String', ...
        '%s does not contain exactly %d nonempty strings.', name, expected_count);
end
end

function ids = decode_text_array(raw, expected_count, name)
if isnumeric(raw)
    raw = char(raw);
end
raw = squeeze(raw);
if ismatrix(raw) && ~isvector(raw)
    if size(raw, 1) == expected_count
        ids = strings(expected_count, 1);
        for index = 1:expected_count
            ids(index) = string(clean_char_text(raw(index, :)));
        end
        return;
    end
    if size(raw, 2) == expected_count
        ids = strings(expected_count, 1);
        for index = 1:expected_count
            ids(index) = string(clean_char_text(raw(:, index).'));
        end
        return;
    end
end

flat = raw(:).';
if expected_count == 1
    ids = string(clean_char_text(flat));
    ids = ids(:);
    return;
end
zero_positions = find(flat == char(0));
pieces = strings(0, 1);
start_position = 1;
for position = [zero_positions, numel(flat) + 1]
    if position > start_position
        item = clean_char_text(flat(start_position:position - 1));
        if ~isempty(item)
            pieces(end + 1, 1) = string(item); %#ok<AGROW>
        end
    end
    start_position = position + 1;
end
if numel(pieces) ~= expected_count
    error('solution:InvalidH5String', ...
        'Unable to decode %d strings from %s.', expected_count, name);
end
ids = pieces;
end

function text_value = decode_text_item(item, name)
if isstring(item) && isscalar(item)
    text_value = item;
elseif ischar(item) || isnumeric(item)
    text_value = string(clean_char_text(char(item(:).')));
else
    error('solution:InvalidH5String', 'Unsupported item in %s.', name);
end
end

function text_value = clean_char_text(text_value)
text_value = char(text_value);
text_value = text_value(:).';
text_value(text_value == char(0)) = [];
text_value = strtrim(text_value);
end

function data = align_predictions(reference, prediction)
n_structures = numel(reference);
if numel(prediction.structure_id) ~= n_structures
    error('solution:StructureSetMismatch', ...
        'Reference and prediction structure counts differ.');
end
reference_ids = string({reference.structure_id}).';
[matched, prediction_index] = ismember(reference_ids, prediction.structure_id);
if ~all(matched) || numel(unique(prediction_index)) ~= n_structures
    error('solution:StructureSetMismatch', ...
        'Reference and prediction structure_id sets differ.');
end

materials = string({reference.material}).';
n_atoms = reshape([reference.n_atoms], [], 1);
energy_ref = reshape([reference.energy_ref_eV], [], 1);
energy_pred = zeros(n_structures, 1);
stress_ref = vertcat(reference.stress_ref_GPa);
stress_pred = zeros(n_structures, 6);
force_ref = cell(n_structures, 1);
force_pred = cell(n_structures, 1);

for ref_index = 1:n_structures
    pred_index = prediction_index(ref_index);
    first_atom = prediction.offsets(pred_index) + 1;
    last_atom = prediction.offsets(pred_index + 1);
    atom_positions = (first_atom:last_atom).';
    predicted_ids = prediction.atom_id(atom_positions);
    if numel(predicted_ids) ~= n_atoms(ref_index) || ...
            numel(unique(predicted_ids)) ~= numel(predicted_ids)
        error('solution:AtomSetMismatch', ...
            'Prediction atom_id set invalid for structure %s.', reference_ids(ref_index));
    end
    [atom_matched, atom_order] = ismember(reference(ref_index).atom_id, predicted_ids);
    if ~all(atom_matched) || numel(unique(atom_order)) ~= n_atoms(ref_index)
        error('solution:AtomSetMismatch', ...
            'Reference and prediction atom_id sets differ for structure %s.', ...
            reference_ids(ref_index));
    end
    energy_pred(ref_index) = prediction.energy_pred_eV(pred_index);
    stress_pred(ref_index, :) = prediction.stress_pred_GPa(pred_index, :);
    force_ref{ref_index} = reference(ref_index).force_ref_eV_A;
    force_pred{ref_index} = prediction.forces_pred_eV_A( ...
        atom_positions(atom_order), :);
end

data = struct();
data.structure_id = reference_ids;
data.material = materials;
data.n_atoms = n_atoms;
data.energy_ref_eV = energy_ref;
data.energy_pred_eV = energy_pred;
data.energy_ref_atom = energy_ref ./ n_atoms;
data.energy_pred_atom = energy_pred ./ n_atoms;
data.stress_ref = stress_ref;
data.stress_pred = stress_pred;
data.force_ref = force_ref;
data.force_pred = force_pred;
end

function groups = compute_all_metrics(data, energy_pred_atom, force_pred, stress_pred)
group_names = {'overall', 'CsSnBr3', 'Cs2SnBr6', 'hetero_N2'};
groups = struct();
for group_index = 1:numel(group_names)
    group_name = group_names{group_index};
    if strcmp(group_name, 'overall')
        selected = true(numel(data.structure_id), 1);
    else
        selected = data.material == string(group_name);
    end
    if ~any(selected)
        error('solution:MissingGroup', 'No structures found for group %s.', group_name);
    end
    force_ref_block = concatenate_cells(data.force_ref(selected));
    force_pred_block = concatenate_cells(force_pred(selected));
    energy_ref_block = data.energy_ref_atom(selected);
    energy_pred_block = energy_pred_atom(selected);
    stress_ref_block = data.stress_ref(selected, :);
    stress_pred_block = stress_pred(selected, :);

    metric = struct();
    metric.count = sum(selected);
    metric.energy_RMSE_meV_atom = root_mean_square( ...
        1000 .* (energy_pred_block - energy_ref_block));
    metric.energy_R2 = r_squared(energy_ref_block, energy_pred_block);
    metric.force_RMSE_eV_A = root_mean_square( ...
        force_pred_block(:) - force_ref_block(:));
    metric.force_R2 = r_squared(force_ref_block(:), force_pred_block(:));
    metric.stress_RMSE_GPa = root_mean_square( ...
        stress_pred_block(:) - stress_ref_block(:));
    metric.stress_R2 = r_squared(stress_ref_block(:), stress_pred_block(:));
    groups.(group_name) = metric;
end
end

function output = append_group_metrics(output, groups)
output.overall = groups.overall;
output.CsSnBr3 = groups.CsSnBr3;
output.Cs2SnBr6 = groups.Cs2SnBr6;
output.hetero_N2 = groups.hetero_N2;
end

function value = root_mean_square(residual)
residual = double(residual(:));
value = sqrt(mean(residual .^ 2));
end

function value = r_squared(reference, prediction)
reference = double(reference(:));
prediction = double(prediction(:));
sst = sum((reference - mean(reference)) .^ 2);
if sst <= 0
    error('solution:UndefinedR2', 'R2 is undefined because SST is zero.');
end
value = 1 - sum((prediction - reference) .^ 2) / sst;
end

function block = concatenate_cells(parts)
if isempty(parts)
    error('solution:EmptyDataBlock', 'Cannot concatenate an empty data block.');
end
block = parts{1};
for index = 2:numel(parts)
    block = [block; parts{index}]; %#ok<AGROW>
end
end

function write_residual_table(path, data)
n_structures = numel(data.structure_id);
energy_error = 1000 .* (data.energy_pred_atom - data.energy_ref_atom);
force_rmse = zeros(n_structures, 1);
stress_rmse = zeros(n_structures, 1);
max_force_error = zeros(n_structures, 1);
for index = 1:n_structures
    force_error = data.force_pred{index} - data.force_ref{index};
    stress_error = data.stress_pred(index, :) - data.stress_ref(index, :);
    force_rmse(index) = root_mean_square(force_error);
    stress_rmse(index) = root_mean_square(stress_error);
    max_force_error(index) = max(abs(force_error(:)));
end

residual_table = table(data.structure_id, data.material, data.n_atoms, ...
    data.energy_ref_eV, data.energy_pred_eV, energy_error, force_rmse, ...
    stress_rmse, max_force_error, 'VariableNames', ...
    {'structure_id', 'material', 'n_atoms', 'energy_ref_eV', ...
    'energy_pred_eV', 'energy_error_meV_atom', 'force_rmse_eV_A', ...
    'stress_rmse_GPa', 'max_force_abs_error_eV_A'});
writetable(residual_table, path);
end

function assessment = assess_training(path, overall_metrics)
history = readtable(path);
required = {'step', 'train_total_loss', 'test_total_loss'};
if ~all(ismember(required, history.Properties.VariableNames))
    error('solution:InvalidTrainingHistory', ...
        'Training history lacks required columns.');
end
steps = double(history.step);
train_loss = double(history.train_total_loss);
test_loss = double(history.test_total_loss);
if any(~isfinite(steps)) || any(~isfinite(train_loss)) || ...
        any(~isfinite(test_loss)) || any(train_loss <= 0) || any(test_loss <= 0)
    error('solution:InvalidTrainingHistory', ...
        'Training history contains invalid values.');
end
[steps, order] = sort(steps);
train_loss = train_loss(order);
test_loss = test_loss(order);
if numel(unique(steps)) ~= numel(steps)
    error('solution:InvalidTrainingHistory', 'Training steps are not unique.');
end

[best_test_loss, best_index] = min(test_loss);
train_reduction = (train_loss(1) - train_loss(end)) / train_loss(1);
test_reduction = (test_loss(1) - test_loss(end)) / test_loss(1);
final_to_best = test_loss(end) / best_test_loss;
final_test_to_train = test_loss(end) / train_loss(end);
converged = train_reduction >= 0.80 && test_reduction >= 0.80;
no_overfitting = final_to_best <= 1.10 && final_test_to_train <= 1.35;
quality_pass = overall_metrics.energy_RMSE_meV_atom <= 2.0 && ...
    overall_metrics.force_RMSE_eV_A <= 0.08 && ...
    overall_metrics.stress_RMSE_GPa <= 0.10 && ...
    overall_metrics.energy_R2 >= 0.99 && ...
    overall_metrics.force_R2 >= 0.70 && ...
    overall_metrics.stress_R2 >= 0.95;

assessment = struct();
assessment.schema_version = '1.0';
assessment.synthetic_data = true;
assessment.n_records = numel(steps);
assessment.first_step = steps(1);
assessment.last_step = steps(end);
assessment.initial_train_total_loss = train_loss(1);
assessment.final_train_total_loss = train_loss(end);
assessment.initial_test_total_loss = test_loss(1);
assessment.final_test_total_loss = test_loss(end);
assessment.train_loss_reduction_fraction = train_reduction;
assessment.test_loss_reduction_fraction = test_reduction;
assessment.best_test_step = steps(best_index);
assessment.final_to_best_test_loss_ratio = final_to_best;
assessment.final_test_to_train_loss_ratio = final_test_to_train;
assessment.converged = converged;
assessment.no_overfitting = no_overfitting;
assessment.metric_quality_pass = quality_pass;
assessment.reliable_for_thermal_transport = ...
    quality_pass && converged && no_overfitting;
end

function [theta, solver_info, convexity] = calibrate_predictions(data)
energy_ref = data.energy_ref_atom(:);
energy_pred = data.energy_pred_atom(:);
force_ref = concatenate_cells(data.force_ref);
force_pred = concatenate_cells(data.force_pred);
stress_ref = data.stress_ref(:);
stress_pred = data.stress_pred(:);
force_ref = force_ref(:);
force_pred = force_pred(:);

n_energy = numel(energy_ref);
n_force = numel(force_ref);
n_stress = numel(stress_ref);
energy_weight = 1 / (0.01 * sqrt(n_energy));
force_weight = 1 / (0.10 * sqrt(n_force));
stress_weight = 1 / (0.10 * sqrt(n_stress));

c_energy = [ones(n_energy, 1), zeros(n_energy, 2)] .* energy_weight;
d_energy = (energy_ref - energy_pred) .* energy_weight;
c_force = [zeros(n_force, 1), force_pred, zeros(n_force, 1)] .* force_weight;
d_force = force_ref .* force_weight;
c_stress = [zeros(n_stress, 2), stress_pred] .* stress_weight;
d_stress = stress_ref .* stress_weight;

lambda = 0.01;
center = [0; 1; 1];
scales = [0.005; 0.05; 0.05];
c_regularization = sqrt(lambda) .* diag(1 ./ scales);
d_regularization = c_regularization * center;
c = [c_energy; c_force; c_stress; c_regularization];
d = [d_energy; d_force; d_stress; d_regularization];

regularization_rank = rank(c_regularization);
augmented_rank = rank(c);
singular_values = svd(c, 'econ');
strictly_convex = regularization_rank == 3 && augmented_rank == 3 && ...
    min(singular_values) > 0;

lower_bounds = [-0.005; 0.80; 0.80];
upper_bounds = [0.005; 1.20; 1.20];
theta0 = center;
options = optimoptions('lsqlin', 'Display', 'off');
[theta, resnorm, ~, exitflag, output] = lsqlin(c, d, [], [], [], [], ...
    lower_bounds, upper_bounds, theta0, options);
if any(~isfinite(theta)) || ~isfinite(resnorm)
    error('solution:CalibrationFailed', 'lsqlin returned non-finite results.');
end

iterations = 0;
if isfield(output, 'iterations')
    iterations = output.iterations;
elseif isfield(output, 'iteration')
    iterations = output.iteration;
end
solver_info = struct('exitflag', exitflag, 'iterations', iterations, ...
    'resnorm', resnorm);
convexity = struct('regularization_rank', regularization_rank, ...
    'augmented_design_rank', augmented_rank, 'n_parameters', 3, ...
    'minimum_singular_value', min(singular_values), ...
    'strictly_convex', strictly_convex);
end

function scaled = scale_cell_matrices(values, scale)
scaled = cell(size(values));
for index = 1:numel(values)
    scaled{index} = scale .* values{index};
end
end

function stats = residual_statistics(data, energy_pred_atom, force_pred, stress_pred)
energy_residual = 1000 .* (energy_pred_atom(:) - data.energy_ref_atom(:));
force_ref = concatenate_cells(data.force_ref);
force_prediction = concatenate_cells(force_pred);
force_residual = force_prediction(:) - force_ref(:);
stress_residual = stress_pred(:) - data.stress_ref(:);

stats = struct();
stats.energy_meV_atom = fit_normal_statistics(energy_residual);
stats.force_eV_A = fit_normal_statistics(force_residual);
stats.stress_GPa = fit_normal_statistics(stress_residual);
end

function result = fit_normal_statistics(residual)
residual = double(residual(:));
if isempty(residual) || any(~isfinite(residual))
    error('solution:InvalidResidual', 'Residual vector is empty or non-finite.');
end
distribution = fitdist(residual, 'Normal');
result = struct();
result.mean = distribution.mu;
result.std = distribution.sigma;
result.rmse = root_mean_square(residual);
result.count = numel(residual);
end

function calibration = build_calibration_document(theta, solver_info, convexity, ...
        metrics_before, metrics_after, stats_before, stats_after)
calibration = struct();
calibration.schema_version = '1.0';
calibration.synthetic_data = true;
calibration.method = 'bounded_regularized_lsqlin';

normalization = struct('energy_eV_atom', 0.01, ...
    'force_eV_A', 0.10, 'stress_GPa', 0.10);
regularization_center = struct('b_E', 0, 's_F', 1, 's_S', 1);
regularization_scales = struct('b_E', 0.005, 's_F', 0.05, 's_S', 0.05);
objective = struct();
objective.normalization_scales = normalization;
objective.energy_normalization_eV_atom = 0.01;
objective.force_normalization_eV_A = 0.10;
objective.stress_normalization_GPa = 0.10;
objective.equal_block_weighting = true;
objective.regularization_lambda = 0.01;
objective.regularization_center = regularization_center;
objective.regularization_scales = regularization_scales;
calibration.objective = objective;

calibration.bounds = struct( ...
    'b_E', struct('lower', -0.005, 'upper', 0.005), ...
    's_F', struct('lower', 0.80, 'upper', 1.20), ...
    's_S', struct('lower', 0.80, 'upper', 1.20));
calibration.parameters = struct( ...
    'energy_bias_eV_atom', theta(1), ...
    'force_scale', theta(2), ...
    'stress_scale', theta(3));
calibration.solver = struct('name', 'lsqlin', ...
    'exitflag', solver_info.exitflag, ...
    'iterations', solver_info.iterations, ...
    'resnorm', solver_info.resnorm);
calibration.objective_value = solver_info.resnorm;
calibration.strict_convexity = convexity;
calibration.unique_solution = convexity.strictly_convex && solver_info.exitflag > 0;
calibration.metrics_before = metrics_before;
calibration.metrics_after = metrics_after;
calibration.residual_statistics_before = stats_before;
calibration.residual_statistics_after = stats_after;
end

function write_json(path, value)
json_text = jsonencode(value);
fid = fopen(path, 'w', 'n', 'UTF-8');
if fid < 0
    error('solution:WriteFailed', 'Unable to write JSON file: %s', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', json_text);
clear cleanup;
end
