function solution(input_dir, output_dir)
input_dir = char(input_dir);
output_dir = char(output_dir);

if ~isfolder(output_dir)
    mkdir(output_dir);
end

systems_path = fullfile(input_dir, 'systems.csv');
trajectory_dir = fullfile(input_dir, 'trajectories');

opts = detectImportOptions(systems_path, 'Delimiter', ',');
text_vars = intersect({'filename', 'material', 'system_type'}, opts.VariableNames);
if ~isempty(text_vars)
    opts = setvartype(opts, text_vars, 'string');
end
systems_table = readtable(systems_path, opts);

n_systems = height(systems_table);
results = repmat(struct( ...
    'material', "", ...
    'filename', "", ...
    'time_ps', [], ...
    'vacf', [], ...
    'frequency_THz', [], ...
    'pdos', [], ...
    'vacf_positive_integral_ps', 0, ...
    'band_fraction_2_3_THz', 0, ...
    'interface_share_2_3_THz', 0, ...
    'dominant_peak_THz', 0, ...
    'band_peak_THz', 0, ...
    'peak_candidates_THz', [], ...
    'interface_band_power', [], ...
    'noninterface_band_power', []), n_systems, 1);

for i = 1:n_systems
    material = string(systems_table.material(i));
    filename = string(systems_table.filename(i));
    file_path = fullfile(trajectory_dir, char(filename));
    results(i) = analyze_system(file_path, material, filename);
end

write_vacf_csv(fullfile(output_dir, 'vacf.csv'), results);
write_pdos_csv(fullfile(output_dir, 'pdos.csv'), results);

interface_stats = build_interface_statistics(results);
comparisons = build_comparisons(results, interface_stats);

if all([ ...
        comparisons.heterostructure_shorter_vacf_than_bulks, ...
        comparisons.heterostructure_enhanced_2_3_THz_fraction, ...
        comparisons.heterostructure_interface_localized_2_3_THz, ...
        comparisons.heterostructure_interface_band_power_significant])
    physical_conclusion = "interface_scattering_shortens_phonon_lifetime_and_localizes_2_3_THz_modes";
else
    physical_conclusion = "";
end

systems_summary = repmat(struct( ...
    'filename', "", ...
    'vacf_positive_integral_ps', 0, ...
    'band_fraction_2_3_THz', 0, ...
    'interface_share_2_3_THz', 0, ...
    'dominant_peak_THz', 0, ...
    'band_peak_THz', 0, ...
    'peak_candidates_THz', []), n_systems, 1);

for i = 1:n_systems
    systems_summary(i).filename = char(results(i).filename);
    systems_summary(i).vacf_positive_integral_ps = finite_scalar(results(i).vacf_positive_integral_ps);
    systems_summary(i).band_fraction_2_3_THz = finite_scalar(results(i).band_fraction_2_3_THz);
    systems_summary(i).interface_share_2_3_THz = finite_scalar(results(i).interface_share_2_3_THz);
    systems_summary(i).dominant_peak_THz = finite_scalar(results(i).dominant_peak_THz);
    systems_summary(i).band_peak_THz = finite_scalar(results(i).band_peak_THz);
    systems_summary(i).peak_candidates_THz = finite_array(results(i).peak_candidates_THz(:).');
end

summary = struct();
summary.schema_version = '1.0';
summary.synthetic_data = true;
summary.systems = systems_summary;
summary.interface_statistics = interface_stats;
summary.comparisons = comparisons;
summary.physical_conclusion = char(physical_conclusion);
summary.method = 'mass-weighted COM removal; time-mean removal; unbiased FFT VACF; symmetric Hann-window unit-integral PDOS; Welch t-test and Hedges g on per-atom 2-3 THz band power';

write_json_utf8(fullfile(output_dir, 'phonon_summary.json'), summary);
end

function result = analyze_system(file_path, material, filename)
time_fs = double(h5read(file_path, '/time_fs'));
time_fs = time_fs(:);
mass_amu = double(h5read(file_path, '/mass_amu'));
mass_amu = mass_amu(:);
region = read_string_vector(h5read(file_path, '/region'));

n_frames = numel(time_fs);
n_atoms = numel(mass_amu);
dt_ps = mean(diff(time_fs)) / 1000;
lag_time_ps = (0:n_frames-1).' * dt_ps;

velocities = read_velocity_array(file_path, n_frames, n_atoms);

mass_weights = reshape(mass_amu, 1, n_atoms, 1);
v_com = sum(velocities .* mass_weights, 2) / sum(mass_amu);
velocities = velocities - v_com;
velocities = velocities - mean(velocities, 1);

signal_matrix = reshape(velocities, n_frames, []);
nfft_vacf = 2 ^ nextpow2(2 * n_frames - 1);
fft_signal = fft(signal_matrix, nfft_vacf, 1);
acf = real(ifft(abs(fft_signal).^2, [], 1));
acf = acf(1:n_frames, :);
acf = acf ./ (n_frames - (0:n_frames-1)).';
vacf = mean(acf, 2);
vacf0 = vacf(1);
if abs(vacf0) > 0
    vacf = vacf / vacf0;
else
    vacf = zeros(n_frames, 1);
    vacf(1) = 1;
end
vacf = finite_array(vacf);
vacf_positive_integral_ps = trapz(lag_time_ps, max(vacf, 0));

window = hann(n_frames, 'symmetric');
windowed_velocities = velocities .* reshape(window, [], 1, 1);
fft_velocity = fft(windowed_velocities, n_frames, 1);
n_freq = floor(n_frames / 2) + 1;
fft_velocity = fft_velocity(1:n_freq, :, :);
frequency_THz = (0:n_freq-1).' / (n_frames * dt_ps);

atomic_power = sum(abs(fft_velocity).^2, 3);
atomic_power = reshape(atomic_power, n_freq, n_atoms);
atomic_power = finite_array(atomic_power);
total_raw_pdos = mean(atomic_power, 2);
total_raw_pdos = finite_array(total_raw_pdos);

pdos_integral = trapz(frequency_THz, total_raw_pdos);
if pdos_integral > 0
    pdos = total_raw_pdos / pdos_integral;
else
    pdos = zeros(size(total_raw_pdos));
end
pdos = finite_array(pdos);

band_fraction_2_3_THz = band_integral(frequency_THz, pdos, 2.0, 3.0);
per_atom_band_power = band_integral(frequency_THz, atomic_power, 2.0, 3.0);
per_atom_band_power = finite_array(per_atom_band_power(:));

interface_mask = strcmpi(region, "interface");
if any(interface_mask)
    total_band_power = sum(per_atom_band_power);
    if total_band_power > 0
        interface_share_2_3_THz = sum(per_atom_band_power(interface_mask)) / total_band_power;
    else
        interface_share_2_3_THz = 0;
    end
else
    interface_share_2_3_THz = 0;
end

positive_mask = frequency_THz > 0.05;
dominant_peak_THz = select_peak_frequency(frequency_THz, total_raw_pdos, positive_mask);
band_peak_THz = select_peak_frequency(frequency_THz, total_raw_pdos, frequency_THz >= 2.0 & frequency_THz <= 3.0);
peak_candidates_THz = top_peak_candidates(frequency_THz, total_raw_pdos, positive_mask, 6);

result = struct();
result.material = material;
result.filename = filename;
result.time_ps = lag_time_ps;
result.vacf = vacf;
result.frequency_THz = frequency_THz;
result.pdos = pdos;
result.vacf_positive_integral_ps = finite_scalar(vacf_positive_integral_ps);
result.band_fraction_2_3_THz = finite_scalar(band_fraction_2_3_THz);
result.interface_share_2_3_THz = finite_scalar(interface_share_2_3_THz);
result.dominant_peak_THz = finite_scalar(dominant_peak_THz);
result.band_peak_THz = finite_scalar(band_peak_THz);
result.peak_candidates_THz = finite_array(peak_candidates_THz(:).');
result.interface_band_power = finite_array(per_atom_band_power(interface_mask));
result.noninterface_band_power = finite_array(per_atom_band_power(~interface_mask));
end

function velocities = read_velocity_array(file_path, n_frames, n_atoms)
raw = double(h5read(file_path, '/velocities_A_ps'));
raw = squeeze(raw);
raw_size = size(raw);
raw_size(end+1:3) = 1;

target_size = [n_frames, n_atoms, 3];
perm_options = perms(1:3);

for i = 1:size(perm_options, 1)
    perm = perm_options(i, :);
    if isequal(raw_size(perm), target_size)
        velocities = permute(raw, perm);
        return;
    end
end

if isequal(raw_size, target_size)
    velocities = raw;
    return;
end

error('Unable to determine velocity dimension order for %s.', file_path);
end

function value = band_integral(frequency, values, low_edge, high_edge)
frequency = frequency(:);
values = double(values);

if high_edge <= low_edge || isempty(frequency) || high_edge < frequency(1) || low_edge > frequency(end)
    if isvector(values)
        value = 0;
    else
        value = zeros(1, size(values, 2));
    end
    return;
end

low_edge = max(low_edge, frequency(1));
high_edge = min(high_edge, frequency(end));
interior = frequency > low_edge & frequency < high_edge;
sample_frequency = [low_edge; frequency(interior); high_edge];
sample_values = interp1(frequency, values, sample_frequency, 'linear');
value = trapz(sample_frequency, sample_values, 1);
value = finite_array(value);
end

function peak_frequency = select_peak_frequency(frequency, power_values, mask)
masked_frequency = frequency(mask);
masked_power = power_values(mask);
if isempty(masked_frequency)
    peak_frequency = 0;
    return;
end
[~, idx] = max(masked_power);
peak_frequency = finite_scalar(masked_frequency(idx));
end

function candidates = top_peak_candidates(frequency, power_values, mask, count)
masked_frequency = frequency(mask);
masked_power = power_values(mask);
if isempty(masked_frequency)
    candidates = zeros(1, 0);
    return;
end

[~, order] = sort(masked_power, 'descend');
count = min(count, numel(order));
candidates = sort(masked_frequency(order(1:count)));
candidates = finite_array(candidates(:).');
end

function interface_stats = build_interface_statistics(results)
materials = string({results.material});
hetero_idx = find(materials == "hetero_N2", 1);

interface_stats = struct();
interface_stats.test = 'Welch two-sample t-test';
interface_stats.alpha = 0.05;
interface_stats.n_interface = 0;
interface_stats.n_noninterface = 0;
interface_stats.mean_interface_band_power = 0;
interface_stats.mean_noninterface_band_power = 0;
interface_stats.mean_difference = 0;
interface_stats.ci95_low = 0;
interface_stats.ci95_high = 0;
interface_stats.t_statistic = 0;
interface_stats.degrees_of_freedom = 0;
interface_stats.p_value = 1;
interface_stats.reject_equal_means = false;
interface_stats.hedges_g = 0;

if isempty(hetero_idx)
    return;
end

interface_power = results(hetero_idx).interface_band_power(:);
noninterface_power = results(hetero_idx).noninterface_band_power(:);

interface_stats.n_interface = numel(interface_power);
interface_stats.n_noninterface = numel(noninterface_power);

if isempty(interface_power) || isempty(noninterface_power)
    return;
end

mean_interface = mean(interface_power);
mean_noninterface = mean(noninterface_power);
mean_difference = mean_interface - mean_noninterface;

[reject_equal_means, p_value, ci95, stats] = ttest2( ...
    interface_power, ...
    noninterface_power, ...
    'Vartype', 'unequal', ...
    'Alpha', 0.05);

n1 = numel(interface_power);
n2 = numel(noninterface_power);
var1 = var(interface_power, 0);
var2 = var(noninterface_power, 0);
pooled_sd = sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / max(n1 + n2 - 2, 1));
if pooled_sd > 0
    cohen_d = mean_difference / pooled_sd;
else
    cohen_d = 0;
end
correction = 1 - 3 / (4 * (n1 + n2) - 9);
hedges_g = correction * cohen_d;

interface_stats.mean_interface_band_power = finite_scalar(mean_interface);
interface_stats.mean_noninterface_band_power = finite_scalar(mean_noninterface);
interface_stats.mean_difference = finite_scalar(mean_difference);
interface_stats.ci95_low = finite_scalar(ci95(1));
interface_stats.ci95_high = finite_scalar(ci95(2));
interface_stats.t_statistic = finite_scalar(stats.tstat);
interface_stats.degrees_of_freedom = finite_scalar(stats.df);
interface_stats.p_value = finite_scalar(p_value);
interface_stats.reject_equal_means = logical(reject_equal_means);
interface_stats.hedges_g = finite_scalar(hedges_g);
end

function comparisons = build_comparisons(results, interface_stats)
materials = string({results.material});
cs_idx = find(materials == "CsSnBr3", 1);
cs2_idx = find(materials == "Cs2SnBr6", 1);
hetero_idx = find(materials == "hetero_N2", 1);

comparisons = struct();
comparisons.heterostructure_shorter_vacf_than_bulks = false;
comparisons.heterostructure_enhanced_2_3_THz_fraction = false;
comparisons.heterostructure_interface_localized_2_3_THz = false;
comparisons.heterostructure_interface_band_power_significant = false;

if isempty(cs_idx) || isempty(cs2_idx) || isempty(hetero_idx)
    return;
end

hetero_vacf = results(hetero_idx).vacf_positive_integral_ps;
bulk_max_fraction = max([results(cs_idx).band_fraction_2_3_THz, results(cs2_idx).band_fraction_2_3_THz]);

comparisons.heterostructure_shorter_vacf_than_bulks = ...
    hetero_vacf < results(cs_idx).vacf_positive_integral_ps && ...
    hetero_vacf < results(cs2_idx).vacf_positive_integral_ps;

comparisons.heterostructure_enhanced_2_3_THz_fraction = ...
    results(hetero_idx).band_fraction_2_3_THz > 2 * bulk_max_fraction;

comparisons.heterostructure_interface_localized_2_3_THz = ...
    results(hetero_idx).interface_share_2_3_THz > 0.90;

comparisons.heterostructure_interface_band_power_significant = ...
    interface_stats.mean_difference > 0 && ...
    interface_stats.ci95_low > 0 && ...
    interface_stats.reject_equal_means;
end

function write_vacf_csv(file_path, results)
fid = fopen(file_path, 'w', 'n', 'UTF-8');
assert(fid >= 0, 'Unable to open %s for writing.', file_path);
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'material,time_ps,vacf\n');
for i = 1:numel(results)
    material = char(results(i).material);
    time_ps = finite_array(results(i).time_ps(:));
    vacf = finite_array(results(i).vacf(:));
    for j = 1:numel(time_ps)
        fprintf(fid, '%s,%.15g,%.15g\n', material, time_ps(j), vacf(j));
    end
end
end

function write_pdos_csv(file_path, results)
fid = fopen(file_path, 'w', 'n', 'UTF-8');
assert(fid >= 0, 'Unable to open %s for writing.', file_path);
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'material,frequency_THz,pdos\n');
for i = 1:numel(results)
    material = char(results(i).material);
    frequency = finite_array(results(i).frequency_THz(:));
    pdos = finite_array(results(i).pdos(:));
    for j = 1:numel(frequency)
        fprintf(fid, '%s,%.15g,%.15g\n', material, frequency(j), pdos(j));
    end
end
end

function write_json_utf8(file_path, value)
json_text = jsonencode(value, 'PrettyPrint', true);
fid = fopen(file_path, 'w', 'n', 'UTF-8');
assert(fid >= 0, 'Unable to open %s for writing.', file_path);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', json_text);
end

function values = finite_array(values)
values = double(values);
values(~isfinite(values)) = 0;
end

function value = finite_scalar(value)
value = double(value);
if ~isfinite(value)
    value = 0;
end
end

function values = read_string_vector(raw)
if isstring(raw)
    values = raw(:);
elseif iscell(raw)
    values = string(raw(:));
elseif ischar(raw)
    if ismatrix(raw) && size(raw, 1) > 1
        values = string(cellstr(raw));
    else
        values = string({raw});
    end
else
    values = string(raw(:));
end
values = strip(values(:));
end
