function check_matlab_ci
%CHECK_MATLAB_CI Exercise the MATLAB capabilities needed by the benchmarks.

root_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(root_dir, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

trace_file = fullfile(result_dir, 'matlab_trace.txt');
diary(trace_file);
trace_cleanup = onCleanup(@() diary('off'));

fprintf('Starting MATLAB capability check.\n');
fprintf('Version: %s\n', version);
fprintf('Release: %s\n', version('-release'));
fprintf('Architecture: %s\n', computer('arch'));

report = struct();
report.generated_at_utc = char(datetime('now', ...
    'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
report.matlab_version = version;
report.release = version('-release');
report.architecture = computer('arch');
report.tests = struct();

checks = {
    'core_numeric', @check_core_numeric;
    'file_io', @() check_file_io(result_dir);
    'signal_processing', @check_signal_processing;
    'statistics', @check_statistics;
    'curve_fitting', @check_curve_fitting;
    'optimization', @check_optimization;
    'parallel_computing', @check_parallel_computing
};

for index = 1:size(checks, 1)
    test_name = checks{index, 1};
    callback = checks{index, 2};
    report.tests.(test_name) = run_check(callback);
    fprintf('%-24s %s\n', test_name, ...
        pass_label(report.tests.(test_name).passed));
end

products = ver;
installed = repmat(struct('name', '', 'version', ''), numel(products), 1);
for index = 1:numel(products)
    installed(index).name = products(index).Name;
    installed(index).version = products(index).Version;
end
report.installed_products = installed;

test_names = fieldnames(report.tests);
passed = false(size(test_names));
for index = 1:numel(test_names)
    passed(index) = report.tests.(test_names{index}).passed;
end
report.passed_count = sum(passed);
report.total_count = numel(passed);
report.all_required_tests_passed = all(passed);

json_file = fullfile(result_dir, 'capability_report.json');
fid = fopen(json_file, 'w');
assert(fid >= 0, 'Cannot create capability_report.json.');
file_cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(report, 'PrettyPrint', true));
clear file_cleanup;

save(fullfile(result_dir, 'capability_report.mat'), 'report');

fprintf('Capability result: %d/%d checks passed.\n', ...
    report.passed_count, report.total_count);
fprintf('All required checks passed: %s\n', ...
    string(report.all_required_tests_passed));
fprintf('Reports written to: %s\n', result_dir);

clear trace_cleanup;
diary('off');

if ~report.all_required_tests_passed
    error('MATLABCapabilityCheck:Failed', ...
        'One or more required MATLAB capability checks failed.');
end
end


function result = run_check(callback)
result = struct('passed', false, 'message', '');
try
    callback();
    result.passed = true;
    result.message = 'ok';
catch exception
    result.message = getReport(exception, 'extended', 'hyperlinks', 'off');
end
end


function check_core_numeric
x = linspace(0, 1, 4096);
y = sin(2 * pi * 17 * x) + 0.2 * cos(2 * pi * 43 * x);
spectrum = fft(y);
assert(numel(spectrum) == numel(y));
assert(all(isfinite(spectrum)));
end


function check_file_io(result_dir)
h5_file = fullfile(result_dir, 'io_test.h5');
csv_file = fullfile(result_dir, 'io_test.csv');
delete_if_present(h5_file);
delete_if_present(csv_file);
cleanup = onCleanup(@() delete_test_files(h5_file, csv_file));

expected = reshape(1:24, [4, 6]);
h5create(h5_file, '/values', size(expected));
h5write(h5_file, '/values', expected);
actual = h5read(h5_file, '/values');
assert(isequal(expected, actual));

table_data = table((1:4)', [0.1; 0.2; 0.3; 0.4], ...
    'VariableNames', {'id', 'value'});
writetable(table_data, csv_file);
loaded = readtable(csv_file);
assert(height(loaded) == 4);

clear cleanup;
delete_test_files(h5_file, csv_file);
end


function check_signal_processing
sample_rate = 200;
time = (0:1/sample_rate:4-1/sample_rate)';
signal = sin(2*pi*12*time) + 0.15*sin(2*pi*37*time);
[power, frequency] = pwelch(signal, [], [], [], sample_rate);
assert(~isempty(power));
assert(numel(power) == numel(frequency));
[correlation, lags] = xcorr(signal, 20, 'coeff');
assert(numel(correlation) == numel(lags));
end


function check_statistics
x = (1:30)';
y = 1.5 + 2.2*x + 0.02*x.^2;
model = fitlm([x, x.^2], y);
assert(model.Rsquared.Ordinary > 0.999);
end


function check_curve_fitting
x = linspace(0, 2, 50)';
y = 0.7 * exp(-1.8*x);
curve = fit(x, y, 'exp1');
predicted = curve(x);
assert(all(isfinite(predicted)));
assert(sqrt(mean((predicted - y).^2)) < 1e-6);
end


function check_optimization
objective = @(x) (x(1)-2)^2 + (x(2)+1)^2;
options = optimoptions('fmincon', 'Display', 'off');
solution = fmincon(objective, [0, 0], [], [], [], [], ...
    [-5, -5], [5, 5], [], options);
assert(norm(solution - [2, -1]) < 1e-4);
end


function check_parallel_computing
existing_pool = gcp('nocreate');
created_pool = isempty(existing_pool);
if created_pool
    pool = parpool('threads', 2);
else
    pool = existing_pool;
end
cleanup = onCleanup(@() close_created_pool(pool, created_pool));
future = parfeval(pool, @sum, 1, 1:100);
value = fetchOutputs(future);
assert(value == 5050);
clear cleanup;
close_created_pool(pool, created_pool);
end


function close_created_pool(pool, created_pool)
if created_pool && ~isempty(pool) && isvalid(pool)
    delete(pool);
end
end


function delete_test_files(h5_file, csv_file)
delete_if_present(h5_file);
delete_if_present(csv_file);
end


function delete_if_present(filename)
if exist(filename, 'file')
    delete(filename);
end
end


function label = pass_label(value)
if value
    label = 'PASS';
else
    label = 'FAIL';
end
end
