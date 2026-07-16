function generate_ground_truth
%GENERATE_GROUND_TRUTH Build deterministic reference outputs for this record.

ground_truth_dir = fileparts(mfilename('fullpath'));
record_dir = fileparts(ground_truth_dir);
input_dir = fullfile(record_dir, 'inputs');
output_dir = fullfile(ground_truth_dir, 'output');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

trace_path = fullfile(output_dir, 'matlab_trace.txt');
if isfile(trace_path)
    delete(trace_path);
end
diary(trace_path);
trace_cleanup = onCleanup(@() diary('off'));

fprintf('NEP/DFT synthetic validation Ground Truth\n');
fprintf('MATLAB version: %s\n', version);
fprintf('MATLAB release: %s\n', version('-release'));
fprintf('Architecture: %s\n', computer('arch'));
fprintf('Synthetic data: true\n');

assert(exist('fitdist', 'file') == 2, ...
    'Statistics and Machine Learning Toolbox is required.');
assert(exist('lsqlin', 'file') == 2, ...
    'Optimization Toolbox is required.');
statistics_toolbox = ver('stats');
optimization_toolbox = ver('optim');
assert(~isempty(statistics_toolbox), ...
    'Cannot query Statistics and Machine Learning Toolbox version.');
assert(~isempty(optimization_toolbox), ...
    'Cannot query Optimization Toolbox version.');
fprintf('Statistics Toolbox version: %s\n', statistics_toolbox.Version);
fprintf('Optimization Toolbox version: %s\n', ...
    optimization_toolbox.Version);

reference_path = fullfile(input_dir, 'test_reference.extxyz');
prediction_path = fullfile(input_dir, 'test_nep_predictions.h5');
history_path = fullfile(input_dir, 'training_history.csv');
assert(isfile(reference_path), 'Missing input: %s', reference_path);
assert(isfile(prediction_path), 'Missing input: %s', prediction_path);
assert(isfile(history_path), 'Missing input: %s', history_path);

fprintf('Reading EXTXYZ reference structures...\n');
frames = readExtxyz(reference_path);
fprintf('Reading HDF5 predictions...\n');
prediction = readPredictions(prediction_path);
fprintf('Reference structures: %d\n', numel(frames));
fprintf('Prediction structures: %d\n', numel(prediction.structure_id));

assert(numel(frames) == 106, 'Expected 106 reference structures.');
reference_ids = {frames.structure_id}.';
assert(numel(unique(reference_ids)) == numel(reference_ids), ...
    'Reference structure_id values must be unique.');
assert(numel(unique(prediction.structure_id)) == ...
    numel(prediction.structure_id), ...
    'Prediction structure_id values must be unique.');
assert(isequal(sort(reference_ids), sort(prediction.structure_id)), ...
    'Reference and prediction structure_id sets differ.');

prediction_index = containers.Map(prediction.structure_id, ...
    num2cell(1:numel(prediction.structure_id)));
n_frames = numel(frames);
aligned_forces = cell(n_frames, 1);
prediction_index_by_frame = zeros(n_frames, 1);

structure_id = strings(n_frames, 1);
material = strings(n_frames, 1);
n_atoms = zeros(n_frames, 1);
energy_ref_eV = zeros(n_frames, 1);
energy_pred_eV = zeros(n_frames, 1);
energy_error_meV_atom = zeros(n_frames, 1);
force_rmse_eV_A = zeros(n_frames, 1);
stress_rmse_GPa = zeros(n_frames, 1);
max_force_abs_error_eV_A = zeros(n_frames, 1);

fprintf('Aligning structures by structure_id and atoms by atom_id...\n');
for frame_index = 1:n_frames
    frame = frames(frame_index);
    pred_index = prediction_index(frame.structure_id);
    prediction_index_by_frame(frame_index) = pred_index;

    first_atom = prediction.offsets(pred_index) + 1;
    last_atom = prediction.offsets(pred_index + 1);
    pred_atom_ids = prediction.atom_id(first_atom:last_atom);
    pred_forces = prediction.forces_eV_A(first_atom:last_atom, :);

    assert(numel(pred_atom_ids) == frame.n_atoms, ...
        'Atom count mismatch for %s.', frame.structure_id);
    assert(numel(unique(pred_atom_ids)) == numel(pred_atom_ids), ...
        'Duplicate prediction atom_id for %s.', frame.structure_id);
    assert(numel(unique(frame.atom_id)) == numel(frame.atom_id), ...
        'Duplicate reference atom_id for %s.', frame.structure_id);

    [matched, row_order] = ismember(frame.atom_id, pred_atom_ids);
    assert(all(matched) && numel(unique(row_order)) == frame.n_atoms, ...
        'atom_id set mismatch for %s.', frame.structure_id);
    aligned = pred_forces(row_order, :);
    aligned_forces{frame_index} = aligned;

    force_error = aligned - frame.forces_ref_eV_A;
    stress_error = prediction.stress_GPa(pred_index, :) - ...
        frame.stress_ref_GPa;

    structure_id(frame_index) = string(frame.structure_id);
    material(frame_index) = string(frame.material);
    n_atoms(frame_index) = frame.n_atoms;
    energy_ref_eV(frame_index) = frame.energy_ref_eV;
    energy_pred_eV(frame_index) = prediction.energy_eV(pred_index);
    energy_error_meV_atom(frame_index) = ...
        (energy_pred_eV(frame_index) - energy_ref_eV(frame_index)) * ...
        1000.0 / frame.n_atoms;
    force_rmse_eV_A(frame_index) = sqrt(mean(force_error(:).^2));
    stress_rmse_GPa(frame_index) = sqrt(mean(stress_error(:).^2));
    max_force_abs_error_eV_A(frame_index) = max(abs(force_error(:)));
end

residuals = table(structure_id, material, n_atoms, energy_ref_eV, ...
    energy_pred_eV, energy_error_meV_atom, force_rmse_eV_A, ...
    stress_rmse_GPa, max_force_abs_error_eV_A);

group_names = {'overall', 'CsSnBr3', 'Cs2SnBr6', 'hetero_N2'};
metrics = struct;
metrics.schema_version = '1.0';
metrics.synthetic_data = true;
metrics.stress_voigt_order = {'xx', 'yy', 'zz', 'yz', 'xz', 'xy'};
for group_index = 1:numel(group_names)
    group_name = group_names{group_index};
    if strcmp(group_name, 'overall')
        mask = true(n_frames, 1);
    else
        mask = strcmp({frames.material}.', group_name);
    end
    metrics.(group_name) = computeGroupMetrics(frames, prediction, ...
        aligned_forces, prediction_index_by_frame, mask);
end

fprintf('Solving bounded regularized posterior calibration with lsqlin...\n');
comparison_data = collectComparisonData(frames, prediction, ...
    aligned_forces, prediction_index_by_frame);
[calibration_parameters, calibration_report] = ...
    solvePosteriorCalibration(comparison_data);

calibrated_metrics = struct;
calibrated_metrics.schema_version = '1.0';
calibrated_metrics.synthetic_data = true;
calibrated_metrics.stress_voigt_order = ...
    {'xx', 'yy', 'zz', 'yz', 'xz', 'xy'};
for group_index = 1:numel(group_names)
    group_name = group_names{group_index};
    if strcmp(group_name, 'overall')
        mask = true(n_frames, 1);
    else
        mask = strcmp({frames.material}.', group_name);
    end
    calibrated_metrics.(group_name) = computeGroupMetrics(frames, ...
        prediction, aligned_forces, prediction_index_by_frame, mask, ...
        calibration_parameters);
end

fprintf('Fitting Normal residual distributions with fitdist...\n');
identity_parameters = struct('energy_bias_eV_atom', 0.0, ...
    'force_scale', 1.0, 'stress_scale', 1.0);
calibration_report.metrics_before = extractMetricGroups(metrics);
calibration_report.metrics_after = extractMetricGroups(calibrated_metrics);
calibration_report.residual_statistics_before = ...
    summarizeResiduals(comparison_data, identity_parameters);
calibration_report.residual_statistics_after = ...
    summarizeResiduals(comparison_data, calibration_parameters);

thresholds = struct;
thresholds.energy_RMSE_meV_atom_max = 2.0;
thresholds.force_RMSE_eV_A_max = 0.08;
thresholds.stress_RMSE_GPa_max = 0.10;
thresholds.energy_R2_min = 0.99;
thresholds.force_R2_min = 0.70;
thresholds.stress_R2_min = 0.95;
thresholds.loss_reduction_fraction_min = 0.80;
thresholds.final_to_best_test_loss_ratio_max = 1.10;
thresholds.final_test_to_train_loss_ratio_max = 1.35;

overall = metrics.overall;
metric_quality_pass = ...
    overall.energy_RMSE_meV_atom <= thresholds.energy_RMSE_meV_atom_max && ...
    overall.force_RMSE_eV_A <= thresholds.force_RMSE_eV_A_max && ...
    overall.stress_RMSE_GPa <= thresholds.stress_RMSE_GPa_max && ...
    overall.energy_R2 >= thresholds.energy_R2_min && ...
    overall.force_R2 >= thresholds.force_R2_min && ...
    overall.stress_R2 >= thresholds.stress_R2_min;

fprintf('Assessing training and test loss histories...\n');
training_assessment = assessTraining(history_path, thresholds, ...
    metric_quality_pass);

metrics_path = fullfile(output_dir, 'metrics.json');
residuals_path = fullfile(output_dir, 'residuals.csv');
training_path = fullfile(output_dir, 'training_assessment.json');
calibration_path = fullfile(output_dir, 'calibration.json');
mat_path = fullfile(output_dir, 'ground_truth.mat');

writeJson(metrics_path, metrics);
writetable(residuals, residuals_path);
writeJson(training_path, training_assessment);
writeJson(calibration_path, calibration_report);

provenance = struct;
provenance.synthetic_data = true;
provenance.generator = 'ground_truth/generate_ground_truth.m';
provenance.input_files = { ...
    'inputs/test_reference.extxyz', ...
    'inputs/test_nep_predictions.h5', ...
    'inputs/training_history.csv'};
provenance.structure_alignment = 'structure_id';
provenance.atom_alignment = 'atom_id';
provenance.statistics_toolbox_function = 'fitdist';
provenance.optimization_toolbox_function = 'lsqlin';
provenance.calibration_method = 'bounded_regularized_lsqlin';
provenance.matlab_release = version('-release');
save(mat_path, 'metrics', 'calibrated_metrics', 'residuals', ...
    'training_assessment', 'calibration_report', ...
    'calibration_parameters', 'thresholds', 'provenance', '-v7');

fprintf('Overall energy RMSE: %.12g meV/atom\n', ...
    overall.energy_RMSE_meV_atom);
fprintf('Overall energy R2: %.12g\n', overall.energy_R2);
fprintf('Overall force RMSE: %.12g eV/Angstrom\n', ...
    overall.force_RMSE_eV_A);
fprintf('Overall force R2: %.12g\n', overall.force_R2);
fprintf('Overall stress RMSE: %.12g GPa\n', ...
    overall.stress_RMSE_GPa);
fprintf('Overall stress R2: %.12g\n', overall.stress_R2);
fprintf('Training converged: %d\n', training_assessment.converged);
fprintf('No overfitting: %d\n', training_assessment.no_overfitting);
fprintf('Reliable for thermal transport: %d\n', ...
    training_assessment.reliable_for_thermal_transport);
fprintf('Calibration b_E: %.12g eV/atom\n', ...
    calibration_parameters.energy_bias_eV_atom);
fprintf('Calibration s_F: %.12g\n', calibration_parameters.force_scale);
fprintf('Calibration s_S: %.12g\n', calibration_parameters.stress_scale);
fprintf('Calibration objective: %.12g\n', ...
    calibration_report.solver.resnorm);
fprintf('Calibration unique solution: %d\n', ...
    calibration_report.unique_solution);
fprintf('Calibrated energy RMSE: %.12g meV/atom\n', ...
    calibrated_metrics.overall.energy_RMSE_meV_atom);
fprintf('Calibrated force RMSE: %.12g eV/Angstrom\n', ...
    calibrated_metrics.overall.force_RMSE_eV_A);
fprintf('Calibrated stress RMSE: %.12g GPa\n', ...
    calibrated_metrics.overall.stress_RMSE_GPa);
fprintf('Wrote Ground Truth outputs to %s\n', output_dir);

diary('off');
clear trace_cleanup;
end

function frames = readExtxyz(path)
fid = fopen(path, 'rt');
assert(fid >= 0, 'Cannot open EXTXYZ file: %s', path);
file_cleanup = onCleanup(@() fclose(fid));

empty_frame = struct( ...
    'structure_id', '', ...
    'material', '', ...
    'n_atoms', 0, ...
    'atom_id', zeros(0, 1), ...
    'forces_ref_eV_A', zeros(0, 3), ...
    'energy_ref_eV', 0, ...
    'stress_ref_GPa', zeros(1, 6));
frames = repmat(empty_frame, 0, 1);

while true
    count_line = fgetl(fid);
    if ~ischar(count_line)
        break;
    end
    count_line = strtrim(count_line);
    if isempty(count_line)
        continue;
    end
    n_atoms = str2double(count_line);
    assert(isfinite(n_atoms) && n_atoms >= 1 && n_atoms == floor(n_atoms), ...
        'Invalid atom count in EXTXYZ.');

    metadata_line = fgetl(fid);
    assert(ischar(metadata_line), 'Missing EXTXYZ metadata line.');
    structure_id = metadataValue(metadata_line, 'structure_id', false);
    material = metadataValue(metadata_line, 'material', false);
    energy_ref = str2double(metadataValue(metadata_line, ...
        'energy_ref_eV', false));
    stress_text = metadataValue(metadata_line, 'stress_ref_GPa', true);
    stress_ref = sscanf(stress_text, '%f').';
    assert(isfinite(energy_ref), 'Invalid reference energy for %s.', ...
        structure_id);
    assert(numel(stress_ref) == 6 && all(isfinite(stress_ref)), ...
        'Invalid reference stress for %s.', structure_id);

    atom_ids = zeros(n_atoms, 1);
    forces = zeros(n_atoms, 3);
    for atom_index = 1:n_atoms
        atom_line = fgetl(fid);
        assert(ischar(atom_line), 'Unexpected end of EXTXYZ frame %s.', ...
            structure_id);
        tokens = strsplit(strtrim(atom_line));
        assert(numel(tokens) == 8, 'Invalid atom row in frame %s.', ...
            structure_id);
        numeric_values = str2double(tokens(2:8));
        assert(all(isfinite(numeric_values)), ...
            'Non-numeric atom row in frame %s.', structure_id);
        atom_ids(atom_index) = numeric_values(4);
        forces(atom_index, :) = numeric_values(5:7);
    end

    frame = empty_frame;
    frame.structure_id = structure_id;
    frame.material = material;
    frame.n_atoms = n_atoms;
    frame.atom_id = atom_ids;
    frame.forces_ref_eV_A = forces;
    frame.energy_ref_eV = energy_ref;
    frame.stress_ref_GPa = stress_ref;
    frames(end + 1, 1) = frame; %#ok<AGROW>
end

clear file_cleanup;
end

function value = metadataValue(line, key, quoted)
escaped_key = regexptranslate('escape', key);
if quoted
    pattern = [escaped_key '="([^"]*)"'];
else
    pattern = [escaped_key '=([^\s]+)'];
end
token = regexp(line, pattern, 'tokens', 'once');
assert(~isempty(token), 'Missing metadata key: %s', key);
value = token{1};
end

function prediction = readPredictions(path)
prediction = struct;
prediction.structure_id = normalizeH5Strings(h5read(path, ...
    '/structure_id'));
n_structures = numel(prediction.structure_id);
prediction.energy_eV = double(h5read(path, '/energy_pred_eV'));
prediction.energy_eV = prediction.energy_eV(:);
prediction.offsets = double(h5read(path, '/offsets'));
prediction.offsets = prediction.offsets(:);
prediction.atom_id = double(h5read(path, '/atom_id'));
prediction.atom_id = prediction.atom_id(:);

assert(numel(prediction.energy_eV) == n_structures, ...
    'HDF5 energy length does not match structure_id length.');
assert(numel(prediction.offsets) == n_structures + 1, ...
    'HDF5 offsets must contain N+1 entries.');
assert(prediction.offsets(1) == 0 && ...
    prediction.offsets(end) == numel(prediction.atom_id) && ...
    all(diff(prediction.offsets) >= 0) && ...
    all(prediction.offsets == floor(prediction.offsets)), ...
    'HDF5 offsets are invalid.');

raw_stress = double(h5read(path, '/stress_pred_GPa'));
prediction.stress_GPa = normalizeComponentMatrix(raw_stress, ...
    n_structures, 6, 'stress_pred_GPa');
raw_forces = double(h5read(path, '/forces_pred_eV_A'));
prediction.forces_eV_A = normalizeComponentMatrix(raw_forces, ...
    numel(prediction.atom_id), 3, 'forces_pred_eV_A');
end

function values = normalizeH5Strings(raw)
if iscell(raw)
    values = cellfun(@decodeH5String, raw(:), 'UniformOutput', false);
elseif isstring(raw)
    values = cellstr(raw(:));
elseif ischar(raw)
    if isrow(raw)
        values = {decodeH5String(raw)};
    else
        values = cellstr(raw);
        values = cellfun(@decodeH5String, values(:), ...
            'UniformOutput', false);
    end
else
    error('Unsupported HDF5 string representation: %s', class(raw));
end
values = values(:);
assert(all(~cellfun(@isempty, values)), ...
    'HDF5 structure_id contains an empty value.');
end

function value = decodeH5String(raw)
if iscell(raw) && isscalar(raw)
    raw = raw{1};
end
if isstring(raw)
    value = char(raw);
elseif ischar(raw)
    value = raw;
elseif isnumeric(raw)
    value = char(raw(:).');
else
    error('Unsupported HDF5 string element: %s', class(raw));
end
value(value == char(0)) = [];
value = strtrim(value);
end

function matrix = normalizeComponentMatrix(raw, n_rows, n_components, name)
raw = squeeze(raw);
if isequal(size(raw), [n_rows, n_components])
    matrix = raw;
elseif isequal(size(raw), [n_components, n_rows])
    matrix = raw.';
else
    error('%s has size %s; expected %d x %d or its transpose.', ...
        name, mat2str(size(raw)), n_rows, n_components);
end
assert(all(isfinite(matrix(:))), '%s contains non-finite values.', name);
end

function result = computeGroupMetrics(frames, prediction, aligned_forces, ...
    prediction_index_by_frame, mask, calibration)
if nargin < 6
    calibration = struct('energy_bias_eV_atom', 0.0, ...
        'force_scale', 1.0, 'stress_scale', 1.0);
end
frame_indices = find(mask);
assert(~isempty(frame_indices), 'Metric group cannot be empty.');
selected_frames = frames(frame_indices);
pred_indices = prediction_index_by_frame(frame_indices);
n_atoms = arrayfun(@(frame) frame.n_atoms, selected_frames);
n_atoms = n_atoms(:);
energy_ref = arrayfun(@(frame) ...
    frame.energy_ref_eV / frame.n_atoms, selected_frames);
energy_ref = energy_ref(:);
energy_pred = prediction.energy_eV(pred_indices) ./ n_atoms + ...
    calibration.energy_bias_eV_atom;
energy_error_meV_atom = (energy_pred - energy_ref) * 1000.0;

force_ref = vertcat(selected_frames.forces_ref_eV_A);
force_pred = vertcat(aligned_forces{frame_indices}) * ...
    calibration.force_scale;
stress_ref = vertcat(selected_frames.stress_ref_GPa);
stress_pred = prediction.stress_GPa(pred_indices, :) * ...
    calibration.stress_scale;

result = struct;
result.count = numel(frame_indices);
result.energy_RMSE_meV_atom = sqrt(mean(energy_error_meV_atom.^2));
result.energy_R2 = canonicalR2(energy_ref, energy_pred);
result.force_RMSE_eV_A = sqrt(mean((force_pred(:) - force_ref(:)).^2));
result.force_R2 = canonicalR2(force_ref(:), force_pred(:));
result.stress_RMSE_GPa = sqrt(mean((stress_pred(:) - stress_ref(:)).^2));
result.stress_R2 = canonicalR2(stress_ref(:), stress_pred(:));
end

function value = canonicalR2(reference, prediction)
reference = double(reference(:));
prediction = double(prediction(:));
residual_sum = sum((prediction - reference).^2);
total_sum = sum((reference - mean(reference)).^2);
if total_sum == 0
    value = double(residual_sum == 0);
    return;
end
value = 1.0 - residual_sum / total_sum;

% Do not replace this definition with squared correlation or a refitted
% regression R-squared; the submitted predictions are not refitted here.
end

function data = collectComparisonData(frames, prediction, aligned_forces, ...
    prediction_index_by_frame)
n_frames = numel(frames);
n_atoms = arrayfun(@(frame) frame.n_atoms, frames);
n_atoms = n_atoms(:);
energy_ref = arrayfun(@(frame) ...
    frame.energy_ref_eV / frame.n_atoms, frames);
energy_ref = energy_ref(:);
energy_pred = prediction.energy_eV(prediction_index_by_frame) ./ n_atoms;

data = struct;
data.energy_ref_eV_atom = energy_ref;
data.energy_pred_eV_atom = energy_pred;
data.force_ref_eV_A = vertcat(frames.forces_ref_eV_A);
data.force_pred_eV_A = vertcat(aligned_forces{1:n_frames});
data.stress_ref_GPa = vertcat(frames.stress_ref_GPa);
data.stress_pred_GPa = prediction.stress_GPa( ...
    prediction_index_by_frame, :);
end

function [parameters, report] = solvePosteriorCalibration(data)
energy_sigma = 0.01;
force_sigma = 0.10;
stress_sigma = 0.10;
regularization_lambda = 0.01;
regularization_center = [0.0; 1.0; 1.0];
regularization_scales = [0.005; 0.05; 0.05];
lower_bounds = [-0.005; 0.80; 0.80];
upper_bounds = [0.005; 1.20; 1.20];

energy_ref = data.energy_ref_eV_atom(:);
energy_pred = data.energy_pred_eV_atom(:);
force_ref = data.force_ref_eV_A(:);
force_pred = data.force_pred_eV_A(:);
stress_ref = data.stress_ref_GPa(:);
stress_pred = data.stress_pred_GPa(:);
n_energy = numel(energy_ref);
n_force = numel(force_ref);
n_stress = numel(stress_ref);

energy_design = [ones(n_energy, 1), zeros(n_energy, 2)] / ...
    (sqrt(n_energy) * energy_sigma);
energy_target = (energy_ref - energy_pred) / ...
    (sqrt(n_energy) * energy_sigma);
force_design = [zeros(n_force, 1), force_pred, ...
    zeros(n_force, 1)] / (sqrt(n_force) * force_sigma);
force_target = force_ref / (sqrt(n_force) * force_sigma);
stress_design = [zeros(n_stress, 2), stress_pred] / ...
    (sqrt(n_stress) * stress_sigma);
stress_target = stress_ref / (sqrt(n_stress) * stress_sigma);

regularization_matrix = sqrt(regularization_lambda) * ...
    diag(1.0 ./ regularization_scales);
regularization_target = regularization_matrix * ...
    regularization_center;
design = [energy_design; force_design; stress_design; ...
    regularization_matrix];
target = [energy_target; force_target; stress_target; ...
    regularization_target];

gram_matrix = design.' * design;
minimum_eigenvalue = min(eig((gram_matrix + gram_matrix.') / 2.0));
full_column_rank = rank(design) == 3;
strictly_convex = full_column_rank && minimum_eigenvalue > 0;
assert(strictly_convex, ...
    'Augmented calibration objective must be strictly convex.');

options = optimoptions('lsqlin', 'Display', 'off');
[theta, resnorm, ~, exitflag, solver_output] = lsqlin( ...
    design, target, [], [], [], [], lower_bounds, upper_bounds, ...
    regularization_center, options);
assert(exitflag > 0, 'lsqlin failed with exitflag %d.', exitflag);

parameters = struct;
parameters.energy_bias_eV_atom = theta(1);
parameters.force_scale = theta(2);
parameters.stress_scale = theta(3);

if isfield(solver_output, 'iterations')
    iterations = solver_output.iterations;
else
    iterations = 0;
end

report = struct;
report.schema_version = '1.0';
report.synthetic_data = true;
report.method = 'bounded_regularized_lsqlin';
report.objective = struct( ...
    'energy_scale_eV_atom', energy_sigma, ...
    'force_scale_eV_A', force_sigma, ...
    'stress_scale_GPa', stress_sigma, ...
    'equal_block_weighting', true, ...
    'regularization_lambda', regularization_lambda, ...
    'regularization_center', regularization_center.', ...
    'regularization_scales', regularization_scales.');
report.bounds = struct( ...
    'energy_bias_eV_atom', [lower_bounds(1), upper_bounds(1)], ...
    'force_scale', [lower_bounds(2), upper_bounds(2)], ...
    'stress_scale', [lower_bounds(3), upper_bounds(3)]);
report.parameters = parameters;
report.solver = struct('name', 'lsqlin', 'exitflag', exitflag, ...
    'iterations', iterations, 'resnorm', resnorm);
report.strict_convexity = struct( ...
    'augmented_design_rank', rank(design), ...
    'parameter_count', 3, ...
    'minimum_gram_eigenvalue', minimum_eigenvalue, ...
    'verified', logical(strictly_convex));
report.unique_solution = logical(strictly_convex && exitflag > 0);
end

function groups = extractMetricGroups(metrics)
groups = struct;
groups.overall = metrics.overall;
groups.CsSnBr3 = metrics.CsSnBr3;
groups.Cs2SnBr6 = metrics.Cs2SnBr6;
groups.hetero_N2 = metrics.hetero_N2;
end

function statistics = summarizeResiduals(data, calibration)
energy_residual = (data.energy_pred_eV_atom + ...
    calibration.energy_bias_eV_atom - data.energy_ref_eV_atom) * 1000.0;
force_residual = calibration.force_scale * data.force_pred_eV_A - ...
    data.force_ref_eV_A;
stress_residual = calibration.stress_scale * data.stress_pred_GPa - ...
    data.stress_ref_GPa;

statistics = struct;
statistics.energy_meV_atom = normalResidualSummary(energy_residual(:));
statistics.force_eV_A = normalResidualSummary(force_residual(:));
statistics.stress_GPa = normalResidualSummary(stress_residual(:));
end

function summary = normalResidualSummary(values)
values = double(values(:));
assert(numel(values) >= 2 && all(isfinite(values)), ...
    'Residual sample is invalid.');
distribution = fitdist(values, 'Normal');
summary = struct;
summary.mean = distribution.mu;
summary.std = distribution.sigma;
summary.rmse = sqrt(mean(values.^2));
summary.count = numel(values);
end

function assessment = assessTraining(path, thresholds, metric_quality_pass)
history = readtable(path, 'VariableNamingRule', 'preserve');
required = {'step', 'train_total_loss', 'test_total_loss'};
assert(all(ismember(required, history.Properties.VariableNames)), ...
    'training_history.csv is missing a required column.');

steps = double(history.step);
train_loss = double(history.train_total_loss);
test_loss = double(history.test_total_loss);
assert(~isempty(steps) && all(isfinite(steps)) && ...
    all(isfinite(train_loss)) && all(isfinite(test_loss)), ...
    'Training history contains invalid values.');

[best_test_loss, best_index] = min(test_loss);
train_reduction = 1.0 - train_loss(end) / train_loss(1);
test_reduction = 1.0 - test_loss(end) / test_loss(1);
final_to_best = test_loss(end) / best_test_loss;
final_test_to_train = test_loss(end) / train_loss(end);

converged = ...
    train_reduction >= thresholds.loss_reduction_fraction_min && ...
    test_reduction >= thresholds.loss_reduction_fraction_min;
no_overfitting = ...
    final_to_best <= thresholds.final_to_best_test_loss_ratio_max && ...
    final_test_to_train <= thresholds.final_test_to_train_loss_ratio_max;

assessment = struct;
assessment.schema_version = '1.0';
assessment.synthetic_data = true;
assessment.n_records = height(history);
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
assessment.converged = logical(converged);
assessment.no_overfitting = logical(no_overfitting);
assessment.metric_quality_pass = logical(metric_quality_pass);
assessment.reliable_for_thermal_transport = logical( ...
    converged && no_overfitting && metric_quality_pass);
assessment.reliability_thresholds = thresholds;
end

function writeJson(path, value)
text = jsonencode(value, 'PrettyPrint', true);
fid = fopen(path, 'wt');
assert(fid >= 0, 'Cannot open JSON output: %s', path);
file_cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', text);
clear file_cleanup;
end
