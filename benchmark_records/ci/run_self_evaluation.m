function run_self_evaluation
%RUN_SELF_EVALUATION Score each Ground Truth output as a reference submission.

project_root = fileparts(fileparts(mfilename('fullpath')));
record_ids = {
    'matlab_gpumd_conductivity_robustness'
    'matlab_nep_dft_energy_force_stress_validation'
    'matlab_vacf_pdos_interface_localization'
};

original_dir = pwd;
cleanup = onCleanup(@() cd(original_dir));
scores = zeros(numel(record_ids), 1);

for index = 1:numel(record_ids)
    record_id = record_ids{index};
    record_dir = fullfile(project_root, 'records', record_id);
    harness_dir = fullfile(record_dir, 'harness');
    output_dir = fullfile(record_dir, 'ground_truth', 'output');
    report_path = fullfile(output_dir, 'self_evaluation.json');

    assert(isfolder(output_dir), 'Missing Ground Truth output: %s', output_dir);
    assert(isfile(fullfile(harness_dir, 'evaluate_submission.m')), ...
        'Missing evaluator for %s', record_id);

    fprintf('\n[%d/%d] Self-evaluating %s\n', ...
        index, numel(record_ids), record_id);
    cd(harness_dir);
    clear evaluate_submission;
    rehash;
    report = evaluate_submission(output_dir, output_dir, report_path);
    assert(isfield(report, 'score'), 'Evaluator did not return score.');
    scores(index) = report.score;
    fprintf('Score: %.6g/100\n', scores(index));
    assert(abs(scores(index) - 100) < 1e-9, ...
        'Reference self-evaluation failed for %s: %.6g', ...
        record_id, scores(index));
end

clear cleanup;
cd(original_dir);
fprintf('\nAll reference submissions scored 100/100.\n');
end
