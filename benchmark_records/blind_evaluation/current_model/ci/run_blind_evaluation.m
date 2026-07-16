function run_blind_evaluation
%RUN_BLIND_EVALUATION Execute and score three no-context model solutions.

blind_root = fileparts(fileparts(mfilename('fullpath')));
benchmark_root = fileparts(fileparts(blind_root));
records_root = fullfile(benchmark_root, 'records');
record_ids = {
    'matlab_gpumd_conductivity_robustness'
    'matlab_nep_dft_energy_force_stress_validation'
    'matlab_vacf_pdos_interface_localization'
};

original_dir = pwd;
original_path = path;
cleanup = onCleanup(@() restore_environment(original_dir, original_path));
summary = struct();
summary.schema_version = '1.0';
summary.evaluation_type = 'no_context_current_model_pretest';
summary.matlab_version = version;
summary.matlab_release = version('-release');
summary.records = struct();

for index = 1:numel(record_ids)
    record_id = record_ids{index};
    blind_record_dir = fullfile(blind_root, record_id);
    submission_dir = fullfile(blind_record_dir, 'submission');
    solution_file = fullfile(submission_dir, 'solution.m');
    input_dir = fullfile(records_root, record_id, 'inputs');
    output_dir = fullfile(submission_dir, 'output');
    result_dir = fullfile(blind_record_dir, 'result');
    execution_trace = fullfile(result_dir, 'matlab_execution_trace.txt');
    score_path = fullfile(result_dir, 'score.json');

    if ~isfolder(result_dir)
        mkdir(result_dir);
    end
    if isfolder(output_dir)
        rmdir(output_dir, 's');
    end
    mkdir(output_dir);

    result = struct();
    result.record_id = record_id;
    result.solution_file = solution_file;
    result.execution_success = false;
    result.evaluation_success = false;
    result.score = 0;
    result.error = '';

    diary(execution_trace);
    fprintf('Blind evaluation record: %s\n', record_id);
    fprintf('MATLAB: %s\n', version);
    fprintf('Input directory: %s\n', input_dir);
    fprintf('Submission directory: %s\n', submission_dir);
    fprintf('Output directory: %s\n', output_dir);

    try
        assert(isfile(solution_file), 'Missing blind solution: %s', solution_file);
        addpath(submission_dir);
        clear solution;
        rehash;
        solution(input_dir, output_dir);
        result.execution_success = true;
        fprintf('MODEL_SOLUTION_EXIT=success\n');
    catch exception
        result.error = getReport(exception, 'extended', 'hyperlinks', 'off');
        fprintf(2, 'MODEL_SOLUTION_EXIT=failure\n%s\n', result.error);
    end
    diary off;

    if result.execution_success
        try
            evaluator_dir = fullfile(records_root, record_id, 'harness');
            ground_truth_dir = fullfile(records_root, record_id, ...
                'ground_truth', 'output');
            cd(evaluator_dir);
            clear evaluate_submission;
            rehash;
            report = evaluate_submission(output_dir, ground_truth_dir, score_path);
            result.evaluation_success = true;
            result.score = report.score;
        catch exception
            result.error = getReport(exception, 'extended', 'hyperlinks', 'off');
            write_failure_score(score_path, record_id, result.error);
        end
    else
        write_failure_score(score_path, record_id, result.error);
    end

    fprintf('Blind score for %s: %.6g/100\n', record_id, result.score);
    summary.records.(record_id) = result;
    cd(original_dir);
    path(original_path);
end

scores = zeros(numel(record_ids), 1);
for index = 1:numel(record_ids)
    scores(index) = summary.records.(record_ids{index}).score;
end
summary.scores = scores;
summary.all_perfect = all(abs(scores - 100) < 1e-9);
summary.mean_score = mean(scores);
summary.minimum_score = min(scores);

summary_path = fullfile(blind_root, 'blind_evaluation_summary.json');
write_json(summary_path, summary);
fprintf('Blind evaluation summary: %s\n', summary_path);
fprintf('All blind solutions perfect: %d\n', summary.all_perfect);

clear cleanup;
restore_environment(original_dir, original_path);
end


function write_failure_score(path_value, record_id, error_text)
report = struct();
report.schema_version = '1.0';
report.record_id = record_id;
report.score = 0;
report.max_score = 100;
report.passed = false;
report.execution_or_evaluation_error = error_text;
write_json(path_value, report);
end


function write_json(path_value, value)
fid = fopen(path_value, 'w');
assert(fid >= 0, 'Cannot write JSON file: %s', path_value);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(value, 'PrettyPrint', true));
clear cleanup;
end


function restore_environment(original_dir, original_path)
cd(original_dir);
path(original_path);
end
