function [] = MNN_ReachTask_AzEl(pNum, epochTypes, lambda)
%% MNN_ReachTask_AzEl - MNN computed per condition (24 locations only, no color)
%
% Same as MNN_ReachTask but uses 24 conditions (location only). Covariance is
% computed per location, averaged across trials from both goColors at that location.
% C_bar = mean across 24 conditions.
%
% Steps:
%   1. Load epoched .set files (from ERP_ReachTask) and behavior.
%   2. Assign condition per trial: location only (1..24). Trials from both goColors
%      at the same location are pooled for that condition.
%   3. For each condition c and time t: cov(trials x channels) -> C_ct(c,t).
%   4. C_c(c) = mean_t C_ct(c,t);  C_bar = mean_c C_c(c).
%   5. C_bar_reg = C_bar + lambda*I;  W = C_bar_reg^{-1/2}.
%   6. Apply W, z-score, windsorize. Save to ERPs_<Epoch>_MNN_AzEl/.
%
% Run after ERP_ReachTask. RDA_ReachTask_AzEl uses these outputs.

if nargin < 1
    pNum = input('Participant Number: ');
end
if nargin < 2 || isempty(epochTypes)
    epochTypes = {'Target','GoCue','Movement','Touch'};
end
if nargin < 3 || isempty(lambda)
    lambda = [];  % set below from data
end

rootPath = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
behavGoodDir = fullfile(rootPath, 'processed', sprintf('p%d', pNum), 'matlab_files_good');

if ~exist(behavGoodDir, 'dir')
    error('matlab_files_good not found: %s (run Preprocess_ReachTask first)', behavGoodDir);
end

%% Process each epoch type
for e = 1:numel(epochTypes)
    epName = epochTypes{e};
    inDir  = fullfile(rootPath, 'processed', sprintf('p%d', pNum), ['ERPs_', epName]);
    outDir = fullfile(rootPath, 'processed', sprintf('p%d', pNum), [sprintf('ERPs_%s_MNN_AzEl', epName)]);

    if ~exist(inDir, 'dir')
        warning('Skipping %s: no folder %s', epName, inDir);
        continue;
    end

    files = dir(fullfile(inDir, ['*_', epName, 'Epoch.set']));
    if isempty(files)
        warning('Skipping %s: no *%sEpoch.set in %s', epName, epName, inDir);
        continue;
    end

    % Sort by session, block
    tok = regexp({files.name}, 'p\d+s(\d+)b(\d+)_', 'tokens');
    valid = ~cellfun(@isempty, tok);
    files = files(valid);
    snums = cellfun(@(x) str2double(x{1}{1}), tok(valid));
    bnums = cellfun(@(x) str2double(x{1}{2}), tok(valid));
    [~, ord] = sortrows([snums(:), bnums(:)], [1 2]);
    files = files(ord);

    fprintf('\n=== MNN AzEl (24 conditions): %s (%d blocks) ===\n', epName, numel(files));

    % ---- Pool data and condition labels (24 conditions = location only) ----
    allData = [];
    allCond = [];

    for fIdx = 1:numel(files)
        fname = files(fIdx).name;
        EEG = pop_loadset(fname, inDir);
        [~, base, ~] = fileparts(fname);
        base = strrep(base, ['_', epName, 'Epoch'], '');
        tokB = regexp(base, 'p\d+s(\d+)b(\d+)', 'tokens');
        if isempty(tokB)
            warning('Skip %s: cannot parse session/block', fname);
            continue;
        end
        behFile = fullfile(behavGoodDir, [base, '.mat']);
        if ~exist(behFile, 'file')
            warning('Skip %s: no behavior %s', fname, behFile);
            continue;
        end
        S = load(behFile);
        if isfield(S, 'theData'), D = S.theData; else, D = S.data; end
        locIdx  = D.locIdx(:);
        nTr = size(EEG.data, 3);
        if numel(locIdx) ~= nTr
            warning('Skip %s: behavior trials (%d) ~= EEG trials (%d)', fname, numel(locIdx), nTr);
            continue;
        end
        % Map experiment loc (1..12, 14..25) to theoretical 1..24
        locTheory = locIdx - (locIdx >= 13);
        condIdx = locTheory;  % 24 conditions: location only
        good = ~(isnan(condIdx) | condIdx < 1 | condIdx > 24);
        allData = cat(3, allData, EEG.data(:, :, good));
        allCond = [allCond; condIdx(good)];
    end

    if isempty(allData)
        warning('No data pooled for %s', epName);
        continue;
    end

    [nChans, nTime, nTrials] = size(allData);
    nC = 24;

    % ---- Covariance per condition (location) and time ----
    C_c = zeros(nChans, nChans, nC);
    for c = 1:nC
        idx = find(allCond == c);
        if numel(idx) < 2
            continue;
        end
        C_ct = zeros(nChans, nChans, nTime);
        for t = 1:nTime
            X = squeeze(allData(:, t, idx))';  % (nTrials_c x channels)
            X = X - mean(X, 1);
            C_ct(:, :, t) = (X' * X) / (size(X, 1) - 1);
        end
        C_c(:, :, c) = mean(C_ct, 3);
    end

    % Average across conditions (24 locations)
    C_bar = mean(C_c, 3);
    lam = lambda;
    if isempty(lam)
        lam = 1e-4 * mean(diag(C_bar));
    end
    C_bar = C_bar + lam * eye(nChans);

    % Whitening matrix W = C_bar^{-1/2}
    [V, D] = eig(C_bar);
    d = diag(D);
    d(d < 1e-10) = 1e-10;
    W = V * diag(1 ./ sqrt(d)) * V';
    W = real(W);

    fprintf('  C_bar from 24 conditions, regularized (lambda=%.2e), W applied.\n', lam);

    % ---- Apply W to each block and save ----
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    for fIdx = 1:numel(files)
        fname = files(fIdx).name;
        EEG = pop_loadset(fname, inDir);
        [~, base, ~] = fileparts(fname);
        base = strrep(base, ['_', epName, 'Epoch'], '');
        tokB = regexp(base, 'p\d+s(\d+)b(\d+)', 'tokens');
        if isempty(tokB)
            continue;
        end
        behFile = fullfile(behavGoodDir, [base, '.mat']);
        if ~exist(behFile, 'file')
            continue;
        end
        % Apply MNN: at each time t, data(:,t,:) = W * data(:,t,:)
        for t = 1:size(EEG.data, 2)
            EEG.data(:, t, :) = W * squeeze(EEG.data(:, t, :));
        end

        % Z-score per channel (mean/std across time and trials within this block)
        for ch = 1:size(EEG.data, 1)
            x = EEG.data(ch, :, :);
            mu = mean(x(:));
            sigma = std(x(:));
            if sigma < 1e-10
                sigma = 1;
            end
            EEG.data(ch, :, :) = (EEG.data(ch, :, :) - mu) / sigma;
        end

        % Windsorize: clip to ±3 SD
        EEG.data(EEG.data > 3)  = 3;
        EEG.data(EEG.data < -3) = -3;

        outName = [base, '_', epName, 'Epoch_MNN_AzEl.set'];
        pop_saveset(EEG, 'filename', outName, 'filepath', outDir);
    end

    fprintf('  Saved %d files to %s\n', numel(files), outDir);
end

fprintf('\nMNN AzEl done. Outputs in ERPs_<Epoch>_MNN_AzEl/\n');
end
