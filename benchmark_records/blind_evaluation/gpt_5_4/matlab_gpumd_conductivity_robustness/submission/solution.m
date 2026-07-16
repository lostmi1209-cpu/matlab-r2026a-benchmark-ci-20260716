function solution(input_dir, output_dir)
if isstring(input_dir)
    input_dir = char(input_dir);
end
if isstring(output_dir)
    output_dir = char(output_dir);
end

burn_in_ps = 600;
expected_materials = {'CsSnBr3', 'Cs2SnBr6', 'hetero_N2'};
size_levels = {'small', 'medium', 'large'};
drive_levels = {'low', 'medium', 'high'};
lengths_nm = [8, 16, 23];

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

trace_path = fullfile(output_dir, 'matlab_trace.txt');
trace_fid = fopen(trace_path, 'w');
assert(trace_fid >= 0, 'Failed to open trace file for writing: %s', trace_path);
cleanup_trace = onCleanup(@() fclose(trace_fid)); %#ok<NASGU>

traceLine(trace_fid, 'MATLAB_TRACE_START');
traceLine(trace_fid, 'INPUT_DIR=%s', input_dir);
traceLine(trace_fid, 'OUTPUT_DIR=%s', output_dir);

runs_path = fullfile(input_dir, 'runs.csv');
runs_table = readtable(runs_path, 'FileType', 'text', 'Delimiter', ',', 'ReadVariableNames', true);
required_run_vars = {'run_id', 'material', 'size_label', 'length_nm', ...
    'drive_label', 'driving_force_A_inv', 'seed'};
assert(all(ismember(required_run_vars, runs_table.Properties.VariableNames)), ...
    'runs.csv is missing one or more required columns.');

n_runs = height(runs_table);
assert(n_runs == 81, 'Expected 81 runs, found %d.', n_runs);

run_ids = cell(n_runs, 1);
materials = cell(n_runs, 1);
size_labels = cell(n_runs, 1);
drive_labels = cell(n_runs, 1);
length_nm_col = zeros(n_runs, 1);
driving_force_col = zeros(n_runs, 1);
seed_col = zeros(n_runs, 1);
burn_in_col = burn_in_ps * ones(n_runs, 1);
n_samples_col = zeros(n_runs, 1);
mean_temperature_col = zeros(n_runs, 1);
kappa_in_col = zeros(n_runs, 1);
kappa_out_col = zeros(n_runs, 1);
kappa_total_col = zeros(n_runs, 1);
closure_error_col = zeros(n_runs, 1);

for i = 1:n_runs
    run_id = textAt(runs_table.run_id, i);
    run_ids{i} = run_id;
    materials{i} = textAt(runs_table.material, i);
    size_labels{i} = textAt(runs_table.size_label, i);
    drive_labels{i} = textAt(runs_table.drive_label, i);
    length_nm_col(i) = numericAt(runs_table.length_nm, i);
    driving_force_col(i) = numericAt(runs_table.driving_force_A_inv, i);
    seed_col(i) = numericAt(runs_table.seed, i);

    trajectory_path = fullfile(input_dir, 'raw', run_id, 'kappa.out');
    assert(exist(trajectory_path, 'file') == 2, 'Missing trajectory file: %s', trajectory_path);
    traj = readtable(trajectory_path, 'FileType', 'text', 'Delimiter', ',', 'ReadVariableNames', true);
    required_traj_vars = {'time_ps', 'temperature_K', 'kappa_in_W_mK', ...
        'kappa_out_W_mK', 'kappa_total_W_mK'};
    assert(all(ismember(required_traj_vars, traj.Properties.VariableNames)), ...
        'Trajectory file is missing one or more required columns: %s', trajectory_path);

    time_ps = traj.time_ps;
    assert(all(diff(time_ps) > 0), 'Time axis is not strictly increasing: %s', trajectory_path);

    keep_mask = time_ps >= burn_in_ps;
    n_kept = sum(keep_mask);
    assert(n_kept >= 100, 'Fewer than 100 post-burn-in samples found: %s', trajectory_path);

    temperature_values = traj.temperature_K(keep_mask);
    kappa_in_values = traj.kappa_in_W_mK(keep_mask);
    kappa_out_values = traj.kappa_out_W_mK(keep_mask);
    kappa_total_values = traj.kappa_total_W_mK(keep_mask);

    n_samples_col(i) = n_kept;
    mean_temperature_col(i) = mean(temperature_values);
    kappa_in_col(i) = mean(kappa_in_values);
    kappa_out_col(i) = mean(kappa_out_values);
    kappa_total_col(i) = mean(kappa_total_values);
    closure_error_col(i) = abs(kappa_in_col(i) + kappa_out_col(i) - kappa_total_col(i));
end

assert(numel(unique(run_ids)) == n_runs, 'run_id values must be unique.');
assert(isequal(sort(unique(materials)), sort(expected_materials')), ...
    'Unexpected material set in runs.csv.');
assert(all(isfinite([length_nm_col; driving_force_col; seed_col; burn_in_col; ...
    n_samples_col; mean_temperature_col; kappa_in_col; kappa_out_col; ...
    kappa_total_col; closure_error_col])), 'Encountered non-finite numeric output.');

run_estimates = table(run_ids, materials, size_labels, length_nm_col, drive_labels, ...
    driving_force_col, seed_col, burn_in_col, n_samples_col, mean_temperature_col, ...
    kappa_in_col, kappa_out_col, kappa_total_col, closure_error_col, ...
    'VariableNames', {'run_id', 'material', 'size_label', 'length_nm', ...
    'drive_label', 'driving_force_A_inv', 'seed', 'burn_in_ps', 'n_samples', ...
    'mean_temperature_K', 'kappa_in_W_mK', 'kappa_out_W_mK', ...
    'kappa_total_W_mK', 'component_closure_error_W_mK'});

run_estimates_path = fullfile(output_dir, 'run_estimates.csv');
writetable(run_estimates, run_estimates_path);

normal_reference_coverage = normcdf(1.96) - normcdf(-1.96);
traceLine(trace_fid, 'SYNTHETIC_DATA=true');
traceLine(trace_fid, 'TRAJECTORIES_PROCESSED=%d', n_runs);
traceLine(trace_fid, 'BURN_IN_PS=%d', burn_in_ps);
traceLine(trace_fid, 'POST_BURN_IN_UNIT=independent_run');
traceLine(trace_fid, 'NORMCDF_COVERAGE=%.17g', normal_reference_coverage);

summary = struct();
summary.schema_version = '1.0';
summary.synthetic_data = true;
summary.method = struct();
summary.method.burn_in_ps = burn_in_ps;
summary.method.run_estimator = 'Arithmetic mean of post-burn-in samples for each independent trajectory.';
summary.method.ensemble = 'Material summaries use 27 independent run-level estimates; standard_error is std(kappa_total_W_mK)/sqrt(27) and ci95 is mean +/- 1.96*standard_error.';
summary.method.normal_reference_coverage = normal_reference_coverage;

materials_struct = struct();
size_robustness_struct = struct();
drive_robustness_struct = struct();
size_fit_struct = struct();
material_total_means = zeros(1, numel(expected_materials));

fit_options = optimoptions('lsqcurvefit', 'Display', 'off');
for m = 1:numel(expected_materials)
    material_name = expected_materials{m};
    material_mask = strcmp(run_estimates.material, material_name);
    material_rows = run_estimates(material_mask, :);
    assert(height(material_rows) == 27, 'Expected 27 runs for %s.', material_name);

    total_values = material_rows.kappa_total_W_mK;
    mean_total = mean(total_values);
    standard_error = std(total_values, 0) / sqrt(height(material_rows));
    ci95_low = mean_total - 1.96 * standard_error;
    ci95_high = mean_total + 1.96 * standard_error;
    material_total_means(m) = mean_total;

    material_entry = struct();
    material_entry.kappa_in_W_mK = mean(material_rows.kappa_in_W_mK);
    material_entry.kappa_out_W_mK = mean(material_rows.kappa_out_W_mK);
    material_entry.kappa_total_W_mK = mean_total;
    material_entry.standard_error_W_mK = standard_error;
    material_entry.ci95_low_W_mK = ci95_low;
    material_entry.ci95_high_W_mK = ci95_high;
    material_entry.n_runs = height(material_rows);
    materials_struct.(material_name) = material_entry;

    size_robustness_struct.(material_name) = buildRobustness(material_rows, 'size_label', size_levels);
    drive_robustness_struct.(material_name) = buildRobustness(material_rows, 'drive_label', drive_levels);

    observed_by_size = zeros(1, numel(size_levels));
    for s = 1:numel(size_levels)
        size_mask = strcmp(material_rows.size_label, size_levels{s});
        size_subset = material_rows(size_mask, :);
        assert(height(size_subset) == 9, 'Expected 9 runs for %s at size %s.', material_name, size_levels{s});
        observed_by_size(s) = mean(size_subset.kappa_total_W_mK);
    end

    k_inf0 = max(observed_by_size) + max(max(observed_by_size) - min(observed_by_size), ...
        0.01 * mean(observed_by_size));
    a0 = max(mean((k_inf0 - observed_by_size) .* lengths_nm), eps);

    model_function = @(p, L) p(1) - p(2) ./ L;
    [fit_params, ~, ~, exitflag] = lsqcurvefit(model_function, [k_inf0, a0], ...
        lengths_nm, observed_by_size, [0, 0], [Inf, Inf], fit_options);
    fitted_by_size = model_function(fit_params, lengths_nm);
    rmse_value = sqrt(mean((fitted_by_size - observed_by_size) .^ 2));
    relative_rmse = rmse_value / mean(observed_by_size);
    assert(all(isfinite([fit_params, fitted_by_size, rmse_value, relative_rmse, exitflag])), ...
        'Encountered non-finite fit output for %s.', material_name);
    converged = (exitflag > 0) && (fit_params(1) >= max(observed_by_size)) && ...
        (fit_params(2) >= 0) && (relative_rmse <= 0.01);

    fit_entry = struct();
    fit_entry.model = 'kappa(L)=k_inf-a/L';
    fit_entry.lengths_nm = lengths_nm;
    fit_entry.observed_group_means_W_mK = observed_by_size;
    fit_entry.fitted_group_means_W_mK = fitted_by_size;
    fit_entry.k_inf_W_mK = fit_params(1);
    fit_entry.a_W_nm_mK = fit_params(2);
    fit_entry.rmse_W_mK = rmse_value;
    fit_entry.relative_rmse = relative_rmse;
    fit_entry.exitflag = exitflag;
    fit_entry.converged = logical(converged);
    size_fit_struct.(material_name) = fit_entry;
end

[~, order_idx] = sort(material_total_means, 'descend');
physical_order = expected_materials(order_idx);
expected_order = {'CsSnBr3', 'Cs2SnBr6', 'hetero_N2'};

summary.materials = materials_struct;
summary.physical_order = physical_order;
summary.size_robustness = size_robustness_struct;
summary.drive_robustness = drive_robustness_struct;
summary.size_convergence_fit = size_fit_struct;
summary.conclusions = struct();
summary.conclusions.size_invariant_within_5_percent = all(structfun( ...
    @(s) s.robust_within_5_percent, size_robustness_struct));
summary.conclusions.drive_invariant_within_5_percent = all(structfun( ...
    @(s) s.robust_within_5_percent, drive_robustness_struct));
summary.conclusions.ordering_matches_expected_group_order = all(strcmp(physical_order, expected_order));

json_path = fullfile(output_dir, 'conductivity_summary.json');
json_text = jsonencode(summary);
json_text = strrep(json_text, sprintf('\\/'), '/');
fid_json = fopen(json_path, 'w');
assert(fid_json >= 0, 'Failed to open JSON file for writing: %s', json_path);
cleanup_json = onCleanup(@() fclose(fid_json)); %#ok<NASGU>
fprintf(fid_json, '%s', json_text);

traceLine(trace_fid, 'LSQCURVEFIT_MATERIALS=%d', numel(expected_materials));
traceLine(trace_fid, 'OUTPUT_JSON=conductivity_summary.json');
traceLine(trace_fid, 'OUTPUT_CSV=run_estimates.csv');
traceLine(trace_fid, 'MATLAB_TRACE_END');
end

function robustness = buildRobustness(material_rows, grouping_var, level_order)
group_means = zeros(1, numel(level_order));
group_means_struct = struct();
for i = 1:numel(level_order)
    mask = strcmp(material_rows.(grouping_var), level_order{i});
    subset = material_rows(mask, :);
    assert(height(subset) == 9, 'Expected 9 runs for grouping level %s.', level_order{i});
    group_mean = mean(subset.kappa_total_W_mK);
    group_means(i) = group_mean;
    group_means_struct.(level_order{i}) = group_mean;
end

center_value = mean(group_means);
if center_value == 0
    max_relative_deviation = 0;
else
    max_relative_deviation = max(abs(group_means - center_value) / center_value);
end

robustness = struct();
robustness.group_means_W_mK = group_means_struct;
robustness.max_relative_deviation = max_relative_deviation;
robustness.robust_within_5_percent = logical(max_relative_deviation <= 0.05);
end

function value = textAt(column_data, index)
if isstring(column_data)
    value = char(column_data(index));
elseif iscell(column_data)
    value = char(column_data{index});
elseif ischar(column_data)
    value = strtrim(column_data(index, :));
elseif iscategorical(column_data)
    value = char(string(column_data(index)));
else
    error('Unsupported text column type.');
end
end

function value = numericAt(column_data, index)
if isnumeric(column_data)
    value = double(column_data(index));
elseif iscell(column_data)
    value = str2double(column_data{index});
elseif isstring(column_data)
    value = str2double(column_data(index));
else
    error('Unsupported numeric column type.');
end
assert(isfinite(value), 'Encountered non-finite numeric input value.');
end

function traceLine(fid, varargin)
fprintf(fid, '%s\n', sprintf(varargin{:}));
end
