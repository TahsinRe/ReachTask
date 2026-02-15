function [] = UpdateBadTrials(pNum)
%% UpdateBadTrials - Run after ConvertBDF, BEFORE data cleaning 
%
% For each block:
%   1. Load EEG (.set from eeglab) and behavioral (.mat)
%   2. Verify alignment: # GoCue events in EEG == # trials in theData
%   3. Apply bad-trial criteria
%   4. Save updated theData
%
% Alignment: EEG trials are defined by GoCue (213) event order. Behavioral
% trial t corresponds to the t-th GoCue event. If counts match, alignment is OK.

if nargin < 1
    pNum = input('Participant Number: ');
end


%% Paths (match ConvertBDF)
cd('C:/Users/tahsi/OneDrive/Desktop/EEG/eeglab2024.0'); eeglab; close;
rootPath  = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
eegPath   = fullfile(rootPath, 'processed', sprintf('p%d', pNum), 'eeglab');
behavDir  = fullfile(rootPath, 'data', sprintf('p%d', pNum), 'matlab_files');
saveDir   = fullfile(rootPath, 'processed', sprintf('p%d', pNum), 'matlab_files');

% Fallback: behavior may be in participant folder directly
if ~exist(behavDir, 'dir')
    behavDir = fullfile(rootPath, sprintf('p%d', pNum));
end

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

if ~exist(eegPath, 'dir')
    error('EEG path not found: %s (run ConvertBDF first)', eegPath);
end
if ~exist(behavDir, 'dir')
    error('Behavior path not found: %s', behavDir);
end

%% Bad-trial criteria (adjust as needed)
GO_CUE_TYPE  = '213';

% Reach accuracy: bad if distTarget > mean + 3*SD (per participant, across blocks)
useAccuracy_mean3SD = true;

% RT: mark bad if > 3 SD from median (per participant, across blocks)
useRT_median3SD = true;

% MT: mark bad if > 3 SD from median (mirrors RT)
useMT_median3SD = true;

%% Find blocks that have both EEG and behavior
eegFiles = dir(fullfile(eegPath, '*.set'));
if isempty(eegFiles)
    error('No .set files in %s', eegPath);
end

% Parse block numbers from EEG filenames (e.g. p2s1b5 -> 5)
tok = regexp({eegFiles.name}, 'p\d+s\d+b(\d+)\.set', 'tokens');
blockNums = [];
for k = 1:numel(tok)
    if ~isempty(tok{k})
        blockNums(end+1) = str2double(tok{k}{1}{1}); %#ok<AGROW>
    end
end
blockNums = unique(blockNums);

%% First pass: collect RTs, MTs, distTarget to compute participant stats (for 3SD criterion)
RT_median_ms  = [];
RT_sd_ms      = [];
MT_median_ms  = [];
MT_sd_ms      = [];
distTarget_mean = [];
distTarget_sd   = [];
if useRT_median3SD || useMT_median3SD || useAccuracy_mean3SD
    allRT = [];
    allMT = [];
    allDistTarget = [];
    for iBlock = 1:numel(blockNums)
        b = blockNums(iBlock);
        patBeh = sprintf('p%ds*b%d.mat', pNum, b);
        candBeh = dir(fullfile(behavDir, patBeh));
        if isempty(candBeh), continue; end
        S = load(fullfile(behavDir, candBeh(1).name));
        if isfield(S, 'theData'), D = S.theData; elseif isfield(S, 'data'), D = S.data; else, continue; end
        if useRT_median3SD && isfield(D, 'RT')
            rt = D.RT(:) * 1000;  % s -> ms
            allRT = [allRT; rt(~isnan(rt))]; %#ok<AGROW>
        end
        if useMT_median3SD && isfield(D, 'MT')
            mt = D.MT(:) * 1000;  % s -> ms
            allMT = [allMT; mt(~isnan(mt))]; %#ok<AGROW>
        end
        if useAccuracy_mean3SD && isfield(D, 'distTarget')
            dt = D.distTarget(:);
            allDistTarget = [allDistTarget; dt(~isnan(dt))]; %#ok<AGROW>
        end
    end
    if useAccuracy_mean3SD
        if numel(allDistTarget) >= 10
            distTarget_mean = mean(allDistTarget);
            distTarget_sd   = std(allDistTarget);
            if distTarget_sd < 0.1, distTarget_sd = 0.1; end
            fprintf('Participant %d: distTarget mean = %.1f px, SD = %.1f px (n = %d trials)\n', pNum, distTarget_mean, distTarget_sd, numel(allDistTarget));
        else
            fprintf('Participant %d: too few distTarget (%d) for mean/SD, skipping accuracy criterion.\n', pNum, numel(allDistTarget));
            useAccuracy_mean3SD = false;
        end
    end
    if useRT_median3SD
        if numel(allRT) >= 10
            RT_median_ms = median(allRT);
            RT_sd_ms     = std(allRT);
            if RT_sd_ms < 1, RT_sd_ms = 1; end
            fprintf('Participant %d: RT median = %.0f ms, SD = %.0f ms (n = %d trials)\n', pNum, RT_median_ms, RT_sd_ms, numel(allRT));
        else
            fprintf('Participant %d: too few RTs (%d) for median/SD, skipping RT criterion.\n', pNum, numel(allRT));
            useRT_median3SD = false;
        end
    end
    if useMT_median3SD
        if numel(allMT) >= 10
            MT_median_ms = median(allMT);
            MT_sd_ms     = std(allMT);
            if MT_sd_ms < 1, MT_sd_ms = 1; end
            fprintf('Participant %d: MT median = %.0f ms, SD = %.0f ms (n = %d trials)\n', pNum, MT_median_ms, MT_sd_ms, numel(allMT));
        else
            fprintf('Participant %d: too few MTs (%d) for median/SD, skipping MT criterion.\n', pNum, numel(allMT));
            useMT_median3SD = false;
        end
    end
end

fprintf('Participant %d: found %d EEG blocks. Checking alignment and updating bad trials.\n', pNum, numel(blockNums));

nUpdated = 0;
nSkipped = 0;
nAlignErr = 0;

for iBlock = 1:numel(blockNums)
    b = blockNums(iBlock);

    % Find EEG file for this block
    pat = sprintf('p%ds*b%d.set', pNum, b);
    cand = dir(fullfile(eegPath, pat));
    if isempty(cand)
        continue;
    end
    eegFile = fullfile(eegPath, cand(1).name);

    % Find behavior file for this block
    patBeh = sprintf('p%ds*b%d.mat', pNum, b);
    candBeh = dir(fullfile(behavDir, patBeh));
    if isempty(candBeh)
        fprintf('  Block %d: no behavior file, skip.\n', b);
        nSkipped = nSkipped + 1;
        continue;
    end
    behFile = fullfile(behavDir, candBeh(1).name);

    %% Load EEG and count GoCue events
    EEG = pop_loadset('filename', cand(1).name, 'filepath', eegPath);
    types = eventTypesToStr(EEG);
    goIdx = find(strcmp(types, GO_CUE_TYPE));
    nEEGtrials = numel(goIdx);

    %% Load behavior
    S = load(behFile);
    if isfield(S, 'theData')
        theData = S.theData;
    elseif isfield(S, 'data')
        theData = S.data;
    else
        fprintf('  Block %d: no theData/data in %s, skip.\n', b, candBeh(1).name);
        nSkipped = nSkipped + 1;
        continue;
    end

    nBehTrials = size(theData.locIdx, 1);
    if isfield(theData, 'nTrials')
        nBehTrials = theData.nTrials;
    end

    %% Check alignment
    if nEEGtrials ~= nBehTrials
        fprintf('  Block %d: ALIGNMENT ERROR - EEG has %d GoCue events, behavior has %d trials. SKIP (no update).\n', ...
            b, nEEGtrials, nBehTrials);
        nAlignErr = nAlignErr + 1;
        continue;
    end

    %% Apply bad-trial criteria
    bad = theData.badTrial(:);  % preserve existing (e.g. no-touch)

    % Reach accuracy: bad if distTarget > mean + 3*SD
    if useAccuracy_mean3SD && ~isempty(distTarget_mean) && isfield(theData, 'distTarget')
        dt = theData.distTarget(:);
        bad(~isnan(dt) & (dt > distTarget_mean + 3*distTarget_sd)) = 1;
    end

    % RT: bad if > 3 SD from median
    if isfield(theData, 'RT') && useRT_median3SD && ~isempty(RT_median_ms)
        rt = theData.RT(:) * 1000;  % s -> ms
        bad(~isnan(rt) & (abs(rt - RT_median_ms) > 3*RT_sd_ms)) = 1;
    end

    % MT: bad if > 3 SD from median
    if isfield(theData, 'MT') && useMT_median3SD && ~isempty(MT_median_ms)
        mt = theData.MT(:) * 1000;  % s -> ms
        bad(~isnan(mt) & (abs(mt - MT_median_ms) > 3*MT_sd_ms)) = 1;
    end

    % Missing condition
    if isfield(theData, 'locIdx') && isfield(theData, 'goColor')
        bad(isnan(theData.locIdx) | isnan(theData.goColor)) = 1;
    end

    theData.badTrial = double(bad(:));

    %% Save to processed/p#/matlab_files/
    saveFile = fullfile(saveDir, candBeh(1).name);
    save(saveFile, 'theData', '-v7.3');
    nBad = sum(bad);
    fprintf('  Block %d: aligned (%d trials), %d bad -> saved %s\n', b, nBehTrials, nBad, saveFile);
    nUpdated = nUpdated + 1;
end

fprintf('\nDone: %d blocks updated, %d skipped, %d alignment errors.\n', nUpdated, nSkipped, nAlignErr);
end

%% =====================================================================
function typesStr = eventTypesToStr(EEG)
typesStr = cell(1, numel(EEG.event));
for i = 1:numel(EEG.event)
    t = EEG.event(i).type;
    if isnumeric(t)
        typesStr{i} = num2str(t);
    elseif isstring(t)
        typesStr{i} = char(t);
    elseif ischar(t)
        typesStr{i} = t;
    else
        typesStr{i} = 'UNKNOWN';
    end
end
end
