function run_all_ground_truth
%RUN_ALL_GROUND_TRUTH Generate MATLAB ground truth for all benchmark records.

project_root = fileparts(fileparts(mfilename('fullpath')));
record_ids = {
    'matlab_gpumd_conductivity_robustness'
    'matlab_nep_dft_energy_force_stress_validation'
    'matlab_vacf_pdos_interface_localization'
};

original_dir = pwd;
cleanup = onCleanup(@() cd(original_dir));

fprintf('MATLAB benchmark Ground Truth generation\n');
fprintf('Version: %s\n', version);
fprintf('Release: %s\n', version('-release'));
fprintf('Architecture: %s\n', computer('arch'));
fprintf('Project root: %s\n', project_root);

for index = 1:numel(record_ids)
    record_id = record_ids{index};
    ground_truth_dir = fullfile(project_root, 'records', record_id, ...
        'ground_truth');
    generator_file = fullfile(ground_truth_dir, 'generate_ground_truth.m');
    assert(isfile(generator_file), 'Missing generator: %s', generator_file);

    fprintf('\n[%d/%d] Starting %s\n', index, numel(record_ids), record_id);
    cd(ground_truth_dir);
    clear generate_ground_truth;
    rehash;
    generate_ground_truth();
    fprintf('[%d/%d] Completed %s\n', index, numel(record_ids), record_id);
end

clear cleanup;
cd(original_dir);
fprintf('\nAll Ground Truth generators completed successfully.\n');
end
