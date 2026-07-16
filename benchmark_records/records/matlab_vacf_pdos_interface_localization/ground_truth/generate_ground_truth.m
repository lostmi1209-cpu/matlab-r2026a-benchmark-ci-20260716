% Generate MATLAB Ground Truth for the synthetic VACF/PDOS signal task.
% The script reads only the public HDF5 inputs and does not import legacy
% Python reference outputs.

scriptPath = mfilename('fullpath');
groundTruthDir = fileparts(scriptPath);
recordRoot = fileparts(groundTruthDir);
inputDir = fullfile(recordRoot, 'inputs');
outputDir = fullfile(groundTruthDir, 'output');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

tracePath = fullfile(outputDir, 'matlab_trace.txt');
if isfile(tracePath)
    delete(tracePath);
end
diary(tracePath);
diaryCleanup = onCleanup(@() diary('off'));

fprintf('VACF/PDOS signal-processing Ground Truth generation\n');
fprintf('Started: %s\n', char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX')));
fprintf('MATLAB: %s (%s)\n', version, version('-release'));
fprintf('Input directory: %s\n', inputDir);
fprintf('Output directory: %s\n', outputDir);

try
    if isempty(which('hann'))
        error('VACFPDOS:MissingToolbox', ...
            'hann is unavailable; Signal Processing Toolbox is required.');
    end
    if isempty(which('ttest2'))
        error('VACFPDOS:MissingStatisticsToolbox', ...
            'ttest2 is unavailable; Statistics and Machine Learning Toolbox is required.');
    end
    signalProduct = ver('signal');
    statisticsProduct = ver('stats');
    fprintf('Signal Processing Toolbox: %s\n', signalProduct.Version);
    fprintf('Statistics and Machine Learning Toolbox: %s\n', ...
        statisticsProduct.Version);

    metadataPath = fullfile(inputDir, 'systems.csv');
    if ~isfile(metadataPath)
        error('VACFPDOS:MissingMetadata', 'Missing input file: %s', metadataPath);
    end
    systemsTable = readtable(metadataPath, 'FileType', 'text', ...
        'Delimiter', ',', 'TextType', 'string');
    requiredMetadata = ["filename", "material", "system_type", "layer_count"];
    assertColumns(systemsTable, requiredMetadata, 'systems.csv');
    if height(systemsTable) ~= 3
        error('VACFPDOS:SystemCount', ...
            'Expected three systems, found %d.', height(systemsTable));
    end

    expectedMaterials = ["CsSnBr3"; "Cs2SnBr6"; "hetero_N2"];
    if ~all(ismember(expectedMaterials, systemsTable.material))
        error('VACFPDOS:Materials', ...
            'systems.csv must contain CsSnBr3, Cs2SnBr6, and hetero_N2.');
    end

    vacfMaterial = strings(0, 1);
    vacfTime = zeros(0, 1);
    vacfValue = zeros(0, 1);
    pdosMaterial = strings(0, 1);
    pdosFrequency = zeros(0, 1);
    pdosValue = zeros(0, 1);
    systemsPayload = struct();
    analysisResults = struct();

    for row = 1:height(systemsTable)
        filename = systemsTable.filename(row);
        material = systemsTable.material(row);
        trajectoryPath = fullfile(inputDir, 'trajectories', char(filename));
        if ~isfile(trajectoryPath)
            error('VACFPDOS:MissingTrajectory', ...
                'Missing trajectory: %s', trajectoryPath);
        end

        timeFs = double(h5read(trajectoryPath, '/time_fs'));
        masses = double(h5read(trajectoryPath, '/mass_amu'));
        timeFs = timeFs(:);
        masses = masses(:);
        if numel(timeFs) < 8 || any(~isfinite(timeFs)) || any(diff(timeFs) <= 0)
            error('VACFPDOS:TimeAxis', ...
                'Invalid or non-monotonic time axis: %s', filename);
        end
        if isempty(masses) || any(~isfinite(masses)) || any(masses <= 0)
            error('VACFPDOS:Masses', 'Invalid masses: %s', filename);
        end

        rawVelocity = double(h5read(trajectoryPath, '/velocities_A_ps'));
        velocity = normalizeVelocityDimensions(rawVelocity, numel(timeFs), numel(masses));
        rawRegion = h5read(trajectoryPath, '/region');
        region = normalizeStringDataset(rawRegion, numel(masses));

        dtFs = median(diff(timeFs));
        if max(abs(diff(timeFs) - dtFs)) > max(1.0e-9, abs(dtFs) * 1.0e-8)
            error('VACFPDOS:IrregularSampling', ...
                'Trajectory is not uniformly sampled: %s', filename);
        end
        dtPs = dtFs / 1000.0;
        result = analyzeTrajectory(velocity, masses, region, dtPs);

        key = char(material);
        systemsPayload.(key) = struct( ...
            'filename', char(filename), ...
            'vacf_positive_integral_ps', result.vacf_positive_integral_ps, ...
            'band_fraction_2_3_THz', result.band_fraction_2_3_THz, ...
            'interface_share_2_3_THz', result.interface_share_2_3_THz, ...
            'dominant_peak_THz', result.dominant_peak_THz, ...
            'band_peak_THz', result.band_peak_THz, ...
            'peak_candidates_THz', result.peak_candidates_THz);
        analysisResults.(key) = result;

        nTime = numel(result.time_ps);
        nFrequency = numel(result.frequency_THz);
        vacfMaterial = [vacfMaterial; repmat(material, nTime, 1)]; %#ok<AGROW>
        vacfTime = [vacfTime; result.time_ps]; %#ok<AGROW>
        vacfValue = [vacfValue; result.vacf]; %#ok<AGROW>
        pdosMaterial = [pdosMaterial; repmat(material, nFrequency, 1)]; %#ok<AGROW>
        pdosFrequency = [pdosFrequency; result.frequency_THz]; %#ok<AGROW>
        pdosValue = [pdosValue; result.pdos]; %#ok<AGROW>

        fprintf(['Processed %s: frames=%d, atoms=%d, dt=%.9g ps, ' ...
            'VACF+=%.9g ps, band=%.9g, interface=%.9g\n'], ...
            material, size(velocity, 1), size(velocity, 2), dtPs, ...
            result.vacf_positive_integral_ps, result.band_fraction_2_3_THz, ...
            result.interface_share_2_3_THz);
    end

    vacfTable = table(vacfMaterial, vacfTime, vacfValue, ...
        'VariableNames', {'material', 'time_ps', 'vacf'});
    pdosTable = table(pdosMaterial, pdosFrequency, pdosValue, ...
        'VariableNames', {'material', 'frequency_THz', 'pdos'});

    hetero = systemsPayload.hetero_N2;
    heteroSignal = analysisResults.hetero_N2;
    interfacePower = heteroSignal.band_atomic_power(heteroSignal.interface_mask);
    nonInterfacePower = heteroSignal.band_atomic_power(~heteroSignal.interface_mask);
    if numel(interfacePower) < 2 || numel(nonInterfacePower) < 2
        error('VACFPDOS:StatisticsSampleSize', ...
            'Welch test requires at least two interface and two non-interface atoms.');
    end
    alpha = 0.05;
    [rejectEqualMeans, pValue, confidenceInterval, testStats] = ttest2( ...
        interfacePower, nonInterfacePower, 'Vartype', 'unequal', 'Alpha', alpha);
    nInterface = numel(interfacePower);
    nNonInterface = numel(nonInterfacePower);
    meanInterface = mean(interfacePower);
    meanNonInterface = mean(nonInterfacePower);
    meanDifference = meanInterface - meanNonInterface;
    pooledStandardDeviation = sqrt(((nInterface - 1) * var(interfacePower, 0) + ...
        (nNonInterface - 1) * var(nonInterfacePower, 0)) / ...
        (nInterface + nNonInterface - 2));
    if ~isfinite(pooledStandardDeviation) || pooledStandardDeviation <= 0
        error('VACFPDOS:EffectSize', ...
            'Pooled standard deviation is not positive and finite.');
    end
    cohenD = meanDifference / pooledStandardDeviation;
    correction = 1.0 - 3.0 / (4.0 * (nInterface + nNonInterface) - 9.0);
    hedgesG = correction * cohenD;
    interfaceStatistics = struct( ...
        'test', 'Welch two-sample t-test on per-atom unnormalized 2-3 THz band power', ...
        'alpha', alpha, ...
        'n_interface', nInterface, ...
        'n_noninterface', nNonInterface, ...
        'mean_interface_band_power', meanInterface, ...
        'mean_noninterface_band_power', meanNonInterface, ...
        'mean_difference', meanDifference, ...
        'ci95_low', confidenceInterval(1), ...
        'ci95_high', confidenceInterval(2), ...
        't_statistic', testStats.tstat, ...
        'degrees_of_freedom', testStats.df, ...
        'p_value', pValue, ...
        'reject_equal_means', logical(rejectEqualMeans), ...
        'hedges_g', hedgesG);

    bulkVacf = [systemsPayload.CsSnBr3.vacf_positive_integral_ps, ...
        systemsPayload.Cs2SnBr6.vacf_positive_integral_ps];
    bulkBand = [systemsPayload.CsSnBr3.band_fraction_2_3_THz, ...
        systemsPayload.Cs2SnBr6.band_fraction_2_3_THz];
    comparisons = struct( ...
        'heterostructure_shorter_vacf_than_bulks', ...
            hetero.vacf_positive_integral_ps < min(bulkVacf), ...
        'heterostructure_enhanced_2_3_THz_fraction', ...
            hetero.band_fraction_2_3_THz > 2.0 * max(bulkBand), ...
        'heterostructure_interface_localized_2_3_THz', ...
            hetero.interface_share_2_3_THz > 0.90, ...
        'heterostructure_interface_band_power_significant', ...
            logical(rejectEqualMeans) && meanDifference > 0 && ...
            confidenceInterval(1) > 0);

    summary = struct();
    summary.schema_version = '1.0';
    summary.synthetic_data = true;
    summary.systems = systemsPayload;
    summary.interface_statistics = interfaceStatistics;
    summary.comparisons = comparisons;
    summary.physical_conclusion = ...
        'interface_scattering_shortens_phonon_lifetime_and_localizes_2_3_THz_modes';
    summary.method = ['mass-weighted COM removal; time-mean removal; ' ...
        'unbiased FFT VACF; symmetric Hann-window unit-integral PDOS; ' ...
        'Welch t-test and Hedges g on per-atom 2-3 THz band power'];

    vacfPath = fullfile(outputDir, 'vacf.csv');
    pdosPath = fullfile(outputDir, 'pdos.csv');
    summaryPath = fullfile(outputDir, 'phonon_summary.json');
    writetable(vacfTable, vacfPath);
    writetable(pdosTable, pdosPath);
    writeJson(summaryPath, summary);
    save(fullfile(outputDir, 'ground_truth.mat'), 'summary', 'vacfTable', ...
        'pdosTable', 'analysisResults', 'systemsTable', '-v7');

    fprintf('Wrote: %s\n', vacfPath);
    fprintf('Wrote: %s\n', pdosPath);
    fprintf('Wrote: %s\n', summaryPath);
    fprintf('Wrote: %s\n', fullfile(outputDir, 'ground_truth.mat'));
    fprintf('Completed successfully: %s\n', char(datetime('now', ...
        'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX')));
catch ME
    fprintf(2, 'GROUND-TRUTH GENERATION FAILED: %s\n', ME.message);
    for stackIndex = 1:numel(ME.stack)
        fprintf(2, '  at %s line %d\n', ...
            ME.stack(stackIndex).name, ME.stack(stackIndex).line);
    end
    clear diaryCleanup;
    rethrow(ME);
end

clear diaryCleanup;

function result = analyzeTrajectory(velocity, masses, region, dtPs)
nFrames = size(velocity, 1);
nAtoms = size(velocity, 2);

massShape = reshape(masses, [1, nAtoms, 1]);
centerOfMassVelocity = sum(velocity .* massShape, 2) / sum(masses);
corrected = velocity - centerOfMassVelocity;
centered = corrected - mean(corrected, 1);

transformed = fft(centered, 2 * nFrames, 1);
autocorrelation = real(ifft(abs(transformed).^2, [], 1));
autocorrelation = autocorrelation(1:nFrames, :, :);
unbiasedDenominator = reshape((nFrames:-1:1)', [nFrames, 1, 1]);
autocorrelation = autocorrelation ./ unbiasedDenominator;
vacf = reshape(mean(autocorrelation, [2, 3]), [], 1);
if ~isfinite(vacf(1)) || vacf(1) <= 0
    error('VACFPDOS:VACFZero', 'VACF lag-zero value is not positive.');
end
vacf = vacf / vacf(1);
timePs = (0:nFrames - 1)' * dtPs;

window = hann(nFrames, 'symmetric');
windowed = centered .* reshape(window, [nFrames, 1, 1]);
fullSpectrum = abs(fft(windowed, [], 1)).^2;
nPositive = floor(nFrames / 2) + 1;
positiveSpectrum = fullSpectrum(1:nPositive, :, :);
atomicPower = reshape(sum(positiveSpectrum, 3), [nPositive, nAtoms]);
pdos = mean(atomicPower, 2);
frequencyTHz = (0:nPositive - 1)' / (nFrames * dtPs);
pdosIntegral = trapz(frequencyTHz, pdos);
if ~isfinite(pdosIntegral) || pdosIntegral <= 0
    error('VACFPDOS:PDOSIntegral', 'PDOS integral is not positive.');
end
pdos = pdos / pdosIntegral;

bandMask = frequencyTHz >= 2.0 & frequencyTHz <= 3.0;
if nnz(bandMask) < 2
    error('VACFPDOS:BandResolution', ...
        'Frequency grid has fewer than two points in 2-3 THz.');
end
bandFraction = trapz(frequencyTHz(bandMask), pdos(bandMask));
bandAtomicPower = trapz(frequencyTHz(bandMask), ...
    atomicPower(bandMask, :), 1);
totalBandPower = sum(bandAtomicPower);
interfaceMask = strcmpi(strip(region), "interface");
if any(interfaceMask) && totalBandPower > 0
    interfaceShare = sum(bandAtomicPower(interfaceMask)) / totalBandPower;
else
    interfaceShare = 0.0;
end

vacfPositiveIntegral = trapz(timePs, max(vacf, 0.0));
nonzeroIndices = find(frequencyTHz > 0.05);
if isempty(nonzeroIndices)
    error('VACFPDOS:FrequencyGrid', 'No positive frequency above 0.05 THz.');
end
[~, dominantLocal] = max(pdos(nonzeroIndices));
dominantIndex = nonzeroIndices(dominantLocal);
bandIndices = find(bandMask);
[~, bandPeakLocal] = max(pdos(bandIndices));
bandPeakIndex = bandIndices(bandPeakLocal);
[~, powerOrder] = sort(pdos(nonzeroIndices), 'descend');
candidateCount = min(6, numel(powerOrder));
candidateIndices = sort(nonzeroIndices(powerOrder(1:candidateCount)));

result = struct( ...
    'time_ps', timePs, ...
    'vacf', vacf, ...
    'frequency_THz', frequencyTHz, ...
    'pdos', pdos, ...
    'vacf_positive_integral_ps', vacfPositiveIntegral, ...
    'band_fraction_2_3_THz', bandFraction, ...
    'interface_share_2_3_THz', interfaceShare, ...
    'dominant_peak_THz', frequencyTHz(dominantIndex), ...
    'band_peak_THz', frequencyTHz(bandPeakIndex), ...
    'peak_candidates_THz', frequencyTHz(candidateIndices)', ...
    'band_atomic_power', bandAtomicPower(:), ...
    'interface_mask', interfaceMask(:));
end

function velocity = normalizeVelocityDimensions(rawVelocity, nFrames, nAtoms)
rawVelocity = squeeze(rawVelocity);
dimensions = size(rawVelocity);
if numel(dimensions) ~= 3 || numel(rawVelocity) ~= nFrames * nAtoms * 3
    error('VACFPDOS:VelocityShape', ...
        'Velocity data cannot be interpreted as [frame, atom, xyz].');
end

orders = [1, 2, 3; 1, 3, 2; 2, 1, 3; 2, 3, 1; 3, 1, 2; 3, 2, 1];
velocity = [];
for orderIndex = 1:size(orders, 1)
    order = orders(orderIndex, :);
    if dimensions(order(1)) == nFrames && ...
            dimensions(order(2)) == nAtoms && dimensions(order(3)) == 3
        velocity = permute(rawVelocity, order);
        break;
    end
end
if isempty(velocity)
    error('VACFPDOS:VelocityShape', ...
        'No HDF5 dimension permutation matches [%d, %d, 3].', nFrames, nAtoms);
end
if any(~isfinite(velocity), 'all')
    error('VACFPDOS:VelocityFinite', 'Velocity data contains NaN or Inf.');
end
end

function labels = normalizeStringDataset(raw, expectedCount)
if isstring(raw)
    labels = raw(:);
elseif iscell(raw)
    labels = strings(numel(raw), 1);
    for index = 1:numel(raw)
        value = raw{index};
        if isnumeric(value)
            value = char(value(:)');
        end
        labels(index) = string(value);
    end
elseif ischar(raw)
    if size(raw, 1) == expectedCount
        labels = string(cellstr(raw));
    elseif size(raw, 2) == expectedCount
        labels = string(cellstr(raw'));
    else
        labels = string(raw(:));
    end
elseif isnumeric(raw)
    labels = decodeNumericStrings(raw, expectedCount);
else
    error('VACFPDOS:RegionType', 'Unsupported HDF5 region string type.');
end

labels = strip(erase(labels(:), char(0)));
if numel(labels) ~= expectedCount
    error('VACFPDOS:RegionCount', ...
        'Expected %d region labels, found %d.', expectedCount, numel(labels));
end
end

function labels = decodeNumericStrings(raw, expectedCount)
if size(raw, 2) == expectedCount
    labels = strings(expectedCount, 1);
    for index = 1:expectedCount
        labels(index) = string(char(raw(:, index)'));
    end
elseif size(raw, 1) == expectedCount
    labels = strings(expectedCount, 1);
    for index = 1:expectedCount
        labels(index) = string(char(raw(index, :)));
    end
else
    error('VACFPDOS:RegionShape', ...
        'Numeric region data cannot be decoded into %d labels.', expectedCount);
end
end

function assertColumns(inputTable, requiredNames, sourceName)
actualNames = string(inputTable.Properties.VariableNames);
missing = requiredNames(~ismember(requiredNames, actualNames));
if ~isempty(missing)
    error('VACFPDOS:MissingColumn', 'Missing columns in %s: %s', ...
        sourceName, strjoin(missing, ', '));
end
end

function writeJson(path, payload)
jsonText = jsonencode(payload, 'PrettyPrint', true);
fileId = fopen(path, 'w', 'n', 'UTF-8');
if fileId < 0
    error('VACFPDOS:WriteJSON', 'Cannot open output file: %s', path);
end
fileCleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, '%s\n', jsonText);
clear fileCleanup;
end
