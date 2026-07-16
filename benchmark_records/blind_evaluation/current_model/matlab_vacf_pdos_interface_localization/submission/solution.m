function solution(input_dir, submission_dir)
%SOLUTION Compute VACF, PDOS, and interface-localization statistics.

if nargin ~= 2
    error('solution:InvalidInput', 'Expected input_dir and submission_dir.');
end
input_dir = char(string(input_dir));
submission_dir = char(string(submission_dir));

if ~isfolder(submission_dir)
    mkdir(submission_dir);
end

mapping_file = fullfile(input_dir, 'systems.csv');
mapping = readtable(mapping_file, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
required_columns = ["filename", "material", "system_type", "layer_count"];
if ~all(ismember(required_columns, string(mapping.Properties.VariableNames)))
    error('solution:InvalidMapping', 'systems.csv is missing required columns.');
end

filenames = strip(mapping.filename);
materials = strip(mapping.material);
required_materials = ["CsSnBr3", "Cs2SnBr6", "hetero_N2"];
if height(mapping) ~= 3 || ~all(ismember(required_materials, materials)) || ...
        numel(unique(materials)) ~= 3
    error('solution:InvalidMapping', 'Expected one row for each required material.');
end

results = cell(height(mapping), 1);
for k = 1:height(mapping)
    trajectory_file = fullfile(input_dir, 'trajectories', char(filenames(k)));
    results{k} = analyze_system(trajectory_file, materials(k), filenames(k));
end

write_vacf_csv(fullfile(submission_dir, 'vacf.csv'), results);
write_pdos_csv(fullfile(submission_dir, 'pdos.csv'), results);

hetero_index = find(materials == "hetero_N2", 1);
cs_index = find(materials == "CsSnBr3", 1);
vacancy_index = find(materials == "Cs2SnBr6", 1);
hetero = results{hetero_index};
cs_bulk = results{cs_index};
vacancy_bulk = results{vacancy_index};

interface_power = hetero.atomic_band_power(hetero.interface_mask);
noninterface_power = hetero.atomic_band_power(~hetero.interface_mask);
interface_power = interface_power(:);
noninterface_power = noninterface_power(:);
alpha = 0.05;
[reject, p_value, confidence_interval, test_stats] = ttest2( ...
    interface_power, noninterface_power, 'Vartype', 'unequal', ...
    'Alpha', alpha);

n_interface = numel(interface_power);
n_noninterface = numel(noninterface_power);
mean_interface = mean(interface_power);
mean_noninterface = mean(noninterface_power);
mean_difference = mean_interface - mean_noninterface;
pooled_variance = ((n_interface - 1) * var(interface_power, 0) + ...
    (n_noninterface - 1) * var(noninterface_power, 0)) / ...
    (n_interface + n_noninterface - 2);
pooled_sd = sqrt(pooled_variance);
if isnan(pooled_sd) || isinf(pooled_sd) || pooled_sd <= 0
    error('solution:DegenerateStatistic', 'Pooled standard deviation is not positive.');
end
cohen_d = mean_difference / pooled_sd;
J = 1 - 3 / (4 * (n_interface + n_noninterface) - 9);
hedges_g = J * cohen_d;

interface_statistics = struct();
interface_statistics.test = 'two-sided Welch two-sample t-test';
interface_statistics.alpha = alpha;
interface_statistics.n_interface = n_interface;
interface_statistics.n_noninterface = n_noninterface;
interface_statistics.mean_interface_band_power = mean_interface;
interface_statistics.mean_noninterface_band_power = mean_noninterface;
interface_statistics.mean_difference = mean_difference;
interface_statistics.ci95_low = confidence_interval(1);
interface_statistics.ci95_high = confidence_interval(2);
interface_statistics.t_statistic = test_stats.tstat;
interface_statistics.degrees_of_freedom = test_stats.df;
interface_statistics.p_value = p_value;
interface_statistics.reject_equal_means = logical(reject);
interface_statistics.hedges_g = hedges_g;
assert_finite(struct2array(rmfield(interface_statistics, ...
    {'test', 'reject_equal_means'})), 'interface statistics');

comparisons = struct();
comparisons.heterostructure_shorter_vacf_than_bulks = logical( ...
    hetero.summary.vacf_positive_integral_ps < min( ...
    cs_bulk.summary.vacf_positive_integral_ps, ...
    vacancy_bulk.summary.vacf_positive_integral_ps));
comparisons.heterostructure_enhanced_2_3_THz_fraction = logical( ...
    hetero.summary.band_fraction_2_3_THz > 2 * max( ...
    cs_bulk.summary.band_fraction_2_3_THz, ...
    vacancy_bulk.summary.band_fraction_2_3_THz));
comparisons.heterostructure_interface_localized_2_3_THz = logical( ...
    hetero.summary.interface_share_2_3_THz > 0.90);
comparisons.heterostructure_interface_band_power_significant = logical( ...
    mean_difference > 0 && confidence_interval(1) > 0 && reject);

all_comparisons_hold = comparisons.heterostructure_shorter_vacf_than_bulks && ...
    comparisons.heterostructure_enhanced_2_3_THz_fraction && ...
    comparisons.heterostructure_interface_localized_2_3_THz && ...
    comparisons.heterostructure_interface_band_power_significant;
if all_comparisons_hold
    physical_conclusion = ...
        'interface_scattering_shortens_phonon_lifetime_and_localizes_2_3_THz_modes';
else
    physical_conclusion = 'not_all_interface_scattering_criteria_are_satisfied';
end

systems_summary = struct();
for k = 1:numel(results)
    systems_summary.(char(results{k}.material)) = results{k}.summary;
end

output = struct();
output.schema_version = '1.0';
output.synthetic_data = true;
output.systems = systems_summary;
output.interface_statistics = interface_statistics;
output.comparisons = comparisons;
output.physical_conclusion = physical_conclusion;
output.method = ['mass-weighted COM removal; time-mean removal; ' ...
    'unbiased FFT VACF; symmetric Hann-window unit-integral PDOS; ' ...
    'Welch t-test and Hedges g on per-atom 2-3 THz band power'];
write_json(fullfile(submission_dir, 'phonon_summary.json'), output);
end

function result = analyze_system(trajectory_file, material, filename)
time_fs = double(h5read(trajectory_file, '/time_fs'));
mass_amu = double(h5read(trajectory_file, '/mass_amu'));
velocity_raw = double(h5read(trajectory_file, '/velocities_A_ps'));
region_raw = h5read(trajectory_file, '/region');

time_fs = time_fs(:);
mass_amu = mass_amu(:);
n_frames = numel(time_fs);
n_atoms = numel(mass_amu);
assert_finite(time_fs, 'time');
assert_finite(mass_amu, 'mass');
if n_frames < 2 || n_atoms < 2 || any(mass_amu <= 0)
    error('solution:InvalidTrajectory', 'Invalid time or mass data in %s.', filename);
end

velocity = orient_velocity(velocity_raw, n_frames, n_atoms);
assert_finite(velocity, 'velocity');
regions = normalize_regions(region_raw, n_atoms);

time_ps = time_fs / 1000;
time_steps = diff(time_ps);
dt_ps = mean(time_steps);
if isnan(dt_ps) || isinf(dt_ps) || dt_ps <= 0 || ...
        max(abs(time_steps - dt_ps)) > 1e-8 * max(1, abs(dt_ps))
    error('solution:InvalidTrajectory', 'Time samples must be uniformly increasing.');
end

mass_shape = reshape(mass_amu, 1, n_atoms, 1);
com_velocity = sum(velocity .* mass_shape, 2) / sum(mass_amu);
velocity = velocity - com_velocity;
velocity = velocity - mean(velocity, 1);

channels = reshape(velocity, n_frames, n_atoms * 3);
nfft = 2 ^ nextpow2(2 * n_frames - 1);
channel_fft = fft(channels, nfft, 1);
autocorrelation = real(ifft(abs(channel_fft) .^ 2, nfft, 1));
unbiased_denominator = (n_frames:-1:1).';
vacf = mean(autocorrelation(1:n_frames, :), 2) ./ unbiased_denominator;
if isnan(vacf(1)) || isinf(vacf(1)) || vacf(1) <= 0
    error('solution:DegenerateVACF', 'VACF zero-lag value is not positive.');
end
vacf = vacf / vacf(1);
lag_time_ps = (0:n_frames - 1).' * dt_ps;

window = hann(n_frames, 'symmetric');
windowed_velocity = velocity .* reshape(window, n_frames, 1, 1);
velocity_fft = fft(windowed_velocity, n_frames, 1);
n_frequency = floor(n_frames / 2) + 1;
one_sided_fft = velocity_fft(1:n_frequency, :, :);
atomic_power = reshape(sum(abs(one_sided_fft) .^ 2, 3), ...
    n_frequency, n_atoms);
frequency_THz = (0:n_frequency - 1).' / (n_frames * dt_ps);
total_power = mean(atomic_power, 2);
power_integral = trapz(frequency_THz, total_power);
if isnan(power_integral) || isinf(power_integral) || power_integral <= 0
    error('solution:DegeneratePDOS', 'PDOS integral is not positive.');
end
pdos = total_power / power_integral;

band_mask = frequency_THz >= 2 & frequency_THz <= 3;
if nnz(band_mask) < 2
    error('solution:InsufficientResolution', 'The 2-3 THz band has too few bins.');
end
atomic_band_power = trapz(frequency_THz(band_mask), ...
    atomic_power(band_mask, :), 1);
band_fraction = trapz(frequency_THz(band_mask), pdos(band_mask));
interface_mask = strcmpi(regions, "interface");
if any(interface_mask)
    total_band_power = sum(atomic_band_power);
    if isnan(total_band_power) || isinf(total_band_power) || total_band_power <= 0
        error('solution:DegenerateBand', 'Total band power is not positive.');
    end
    interface_share = sum(atomic_band_power(interface_mask)) / total_band_power;
else
    interface_share = 0;
end

positive_frequency_mask = frequency_THz > 0.05;
positive_indices = find(positive_frequency_mask);
if numel(positive_indices) < 6
    error('solution:InsufficientResolution', 'Fewer than six positive frequency bins.');
end
[~, dominant_offset] = max(pdos(positive_frequency_mask));
dominant_peak = frequency_THz(positive_indices(dominant_offset));
band_indices = find(band_mask);
[~, band_offset] = max(pdos(band_mask));
band_peak = frequency_THz(band_indices(band_offset));
[~, ranked_offsets] = sort(pdos(positive_frequency_mask), 'descend');
candidate_indices = positive_indices(ranked_offsets(1:6));
peak_candidates = sort(frequency_THz(candidate_indices)).';

summary = struct();
summary.filename = char(filename);
summary.vacf_positive_integral_ps = trapz(lag_time_ps, max(vacf, 0));
summary.band_fraction_2_3_THz = band_fraction;
summary.interface_share_2_3_THz = interface_share;
summary.dominant_peak_THz = dominant_peak;
summary.band_peak_THz = band_peak;
summary.peak_candidates_THz = peak_candidates;
assert_finite([summary.vacf_positive_integral_ps, ...
    summary.band_fraction_2_3_THz, summary.interface_share_2_3_THz, ...
    summary.dominant_peak_THz, summary.band_peak_THz, peak_candidates], ...
    'system summary');
assert_finite(vacf, 'VACF');
assert_finite(pdos, 'PDOS');
assert_finite(atomic_band_power, 'atomic band power');

result = struct();
result.material = material;
result.lag_time_ps = lag_time_ps;
result.vacf = vacf;
result.frequency_THz = frequency_THz;
result.pdos = pdos;
result.atomic_band_power = atomic_band_power;
result.interface_mask = interface_mask(:).';
result.summary = summary;
end

function velocity = orient_velocity(raw, n_frames, n_atoms)
raw = squeeze(raw);
raw_size = size(raw);
if numel(raw_size) > 3
    error('solution:InvalidVelocityShape', 'Velocity dataset has too many dimensions.');
end
raw_size(end + 1:3) = 1;
raw = reshape(raw, raw_size);
target_size = [n_frames, n_atoms, 3];
permutations = perms(1:3);
for k = 1:size(permutations, 1)
    order = permutations(k, :);
    if isequal(raw_size(order), target_size)
        velocity = permute(raw, order);
        return;
    end
end
error('solution:InvalidVelocityShape', ...
    'Cannot orient velocity dimensions to [frame, atom, xyz].');
end

function labels = normalize_regions(raw, n_atoms)
if isstring(raw)
    labels = raw(:);
elseif iscell(raw)
    labels = strings(numel(raw), 1);
    for k = 1:numel(raw)
        item = raw{k};
        if ischar(item) || isstring(item)
            labels(k) = string(item);
        elseif isnumeric(item)
            labels(k) = string(char(item(:).'));
        else
            error('solution:InvalidRegion', 'Unsupported cell value in region dataset.');
        end
    end
elseif ischar(raw)
    labels = char_matrix_to_strings(raw, n_atoms);
elseif isnumeric(raw)
    labels = char_matrix_to_strings(char(raw), n_atoms);
else
    error('solution:InvalidRegion', 'Unsupported region dataset representation.');
end

if numel(labels) ~= n_atoms
    error('solution:InvalidRegion', 'Region label count does not match atom count.');
end
for k = 1:numel(labels)
    text = char(labels(k));
    text(text == 0) = [];
    labels(k) = string(strtrim(text));
end
if any(strlength(labels) == 0)
    error('solution:InvalidRegion', 'Empty region label encountered.');
end
labels = labels(:);
end

function labels = char_matrix_to_strings(raw, n_atoms)
if isvector(raw) && n_atoms == 1
    labels = string(raw(:).');
elseif size(raw, 1) == n_atoms
    labels = string(cellstr(raw));
elseif size(raw, 2) == n_atoms
    labels = string(cellstr(raw.'));
else
    error('solution:InvalidRegion', 'Cannot orient region character data.');
end
labels = labels(:);
end

function write_vacf_csv(output_file, results)
fid = open_utf8(output_file);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'material,time_ps,vacf\n');
for k = 1:numel(results)
    item = results{k};
    material = char(item.material);
    for row = 1:numel(item.lag_time_ps)
        fprintf(fid, '%s,%.17g,%.17g\n', material, ...
            item.lag_time_ps(row), item.vacf(row));
    end
end
clear cleanup;
end

function write_pdos_csv(output_file, results)
fid = open_utf8(output_file);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'material,frequency_THz,pdos\n');
for k = 1:numel(results)
    item = results{k};
    material = char(item.material);
    for row = 1:numel(item.frequency_THz)
        fprintf(fid, '%s,%.17g,%.17g\n', material, ...
            item.frequency_THz(row), item.pdos(row));
    end
end
clear cleanup;
end

function write_json(output_file, value)
encoded = jsonencode(value);
fid = open_utf8(output_file);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', encoded);
clear cleanup;
end

function fid = open_utf8(filename)
[fid, message] = fopen(filename, 'w', 'n', 'UTF-8');
if fid < 0
    error('solution:OutputFailure', 'Cannot open %s: %s', filename, message);
end
end

function assert_finite(value, label)
if any(isnan(value(:))) || any(isinf(value(:)))
    error('solution:NonFinite', '%s contains NaN or Inf.', label);
end
end
