function [] = Preprocess_ReachTask(pNum)
%% Preprocess_ReachTask - ICA from session 1, applied block-by-block
%
% Bad trials are dropped from EEG before processing; alignment with behavior is preserved.
%
% 1. Load behavior (processed/p#/matlab_files), drop bad-trial segments from each block
% 2. Concatenate session 1 (good trials only) -> run ICA
% 3. Apply ICA to each block (good trials only), ICLabel, remove bad ICs
% 4. Save: post_ICA/ (EEG), matlab_files_good/ (good-trials-only .mat for 1:1 alignment)
%
% Run after ConvertBDF and UpdateBadTrials.
% Downstream (e.g. Stage4): use post_ICA + matlab_files_good for alignment.

if nargin < 1
    pNum = input('Participant Number: ');
end

%% --- Open EEGLAB ---
cd('C:/Users/tahsi/OneDrive/Desktop/EEG/eeglab2024.0'); eeglab; close;

%% Paths
rootPath   = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
eeglabDir  = 'C:/Users/tahsi/OneDrive/Desktop/EEG/eeglab2024.0';
inDir      = fullfile(rootPath, 'processed', sprintf('p%d', pNum), 'eeglab');
outDir     = fullfile(rootPath, 'processed', sprintf('p%d', pNum), 'post_ICA');
behavDir   = fullfile(rootPath, 'processed', sprintf('p%d', pNum), 'matlab_files');
if ~exist(behavDir, 'dir')
    behavDir = fullfile(rootPath, 'data', sprintf('p%d', pNum), 'matlab_files');
end
behavGoodDir = fullfile(rootPath, 'processed', sprintf('p%d', pNum), 'matlab_files_good');
GO_CUE_TYPE = '213'; 
TRIAL_START_TRIG = '210';   % trial start
TRIAL_END_TRIG   = '216';   % trial end (moveStop)
chanLookup = fullfile(eeglabDir, 'plugins', 'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc');
if ~exist(chanLookup, 'file')
    chanLookup = fullfile(eeglabDir, 'sample_data', 'standard-10-5-cap385.elp');
end

resamprate = 256;
bandpass   = [0.1 40];
nChansUse  = 64;

if ~exist(inDir, 'dir')
    error('Input not found: %s (run ConvertBDF first)', inDir);
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Find blocks (session 1 for ICA, all for application)
files = dir(fullfile(inDir, sprintf('p%ds1b*.set', pNum)));
if isempty(files)
    error('No session 1 blocks in %s', inDir);
end
% Sort by block number
tok = regexp({files.name}, 'b(\d+)\.set', 'tokens');
bnum = cellfun(@(x) str2double(x{1}{1}), tok);
[~, ord] = sort(bnum);
files_s1 = files(ord);

allFiles = dir(fullfile(inDir, '*.set'));
tokA = regexp({allFiles.name}, 's(\d+)b(\d+)\.set', 'tokens');
valid = ~cellfun(@isempty, tokA);
allFiles = allFiles(valid);
tokA = tokA(valid);
snums = cellfun(@(x) str2double(x{1}{1}), tokA);
bnums = cellfun(@(x) str2double(x{1}{2}), tokA);
[~, ordA] = sortrows([snums(:), bnums(:)], [1 2]);
allFiles = allFiles(ordA);

%% ========== PHASE 1: ICA from session 1 (good trials only) ==========
fprintf('Participant %d: Loading session 1 blocks for ICA (dropping bad trials)...\n', pNum);
EEG_all = [];
for k = 1:numel(files_s1)
    EEG = pop_loadset(files_s1(k).name, inDir);
    % Load behavior and drop bad trials from EEG
    tokB = regexp(files_s1(k).name, 's(\d+)b(\d+)\.set', 'tokens');
    if ~isempty(tokB) && exist(behavDir, 'dir')
        sb = dir(fullfile(behavDir, sprintf('p%ds%sb%s.mat', pNum, tokB{1}{1}, tokB{1}{2})));
        if ~isempty(sb)
            S = load(fullfile(behavDir, sb(1).name));
            D = S.theData; if ~isfield(S, 'theData') && isfield(S, 'data'), D = S.data; end
            if isfield(D, 'badTrial')
                [EEG, ~] = dropBadTrialsAndGetGoodBehavior(EEG, D, TRIAL_START_TRIG, TRIAL_END_TRIG);
            end
        end
    end
    EEG = pop_resample(EEG, resamprate);
    EEG = pop_chanedit(EEG, 'lookup', chanLookup);
    EEG = pop_select(EEG, 'channel', 1:nChansUse);
    EEG = pop_basicfilter(EEG, 1:EEG.nbchan, 'Boundary', 'boundary', 'Cutoff', bandpass, ...
        'Design', 'butter', 'Filter', 'bandpass', 'Order', 2, 'RemoveDC', 'on');
    if k == 1
        EEG_all = EEG;
    else
        EEG_all = pop_mergeset(EEG_all, EEG, 1);
    end
end
EEG_all = eeg_checkset(EEG_all);

% Bad channel rejection on concatenated session 1
chanlocs_full = EEG_all.chanlocs;
nChansBefore  = EEG_all.nbchan;
EEG_all = clean_rawdata(EEG_all, 'off', 'off', 0.8, 'off', 'off', 'off');
nChansAfter = EEG_all.nbchan;
% Indices of rejected channels (in original 1:nChansBefore)
labelsKept = {EEG_all.chanlocs.labels};
labelsAll  = {chanlocs_full.labels};
[~, badChanIdx] = setdiff(labelsAll, labelsKept);
badChans_s1 = badChanIdx(:)';
Info.chanlocs = chanlocs_full;
Info.nChansRejected = nChansBefore - nChansAfter;
fprintf('  Session 1: %d channels rejected\n', Info.nChansRejected);

EEG_all = pop_reref(EEG_all, []);

% ICA (PICARD)
fprintf('  Running PICARD ICA on session 1 (%d channels, %d samples)...\n', EEG_all.nbchan, EEG_all.pnts);
nRanks = EEG_all.nbchan;
k_rankdeff = nRanks - rank(EEG_all.data(:, resamprate:resamprate*2));
EEG_all = pop_runica(EEG_all, 'icatype', 'picard', 'pca', max(1, nRanks - k_rankdeff));

% Store ICA solution for applying to other blocks
ICA_weights = EEG_all.icaweights;
ICA_sphere  = EEG_all.icasphere;
clear EEG_all;

%% ========== PHASE 2: Apply ICA block-by-block ==========
fprintf('Participant %d: Applying ICA to %d blocks...\n', pNum, numel(allFiles));

if ~exist(behavGoodDir, 'dir')
    mkdir(behavGoodDir);
end

for k = 1:numel(allFiles)
    fname = allFiles(k).name;
    [~, basename, ~] = fileparts(fname);
    fprintf('  Block %s\n', basename);
    theDataGood = [];

    EEG = pop_loadset(fname, inDir);
    % Load behavior and drop bad trials from EEG
    tokB = regexp(fname, 's(\d+)b(\d+)\.set', 'tokens');
    if ~isempty(tokB) && exist(behavDir, 'dir')
        sb = dir(fullfile(behavDir, sprintf('p%ds%sb%s.mat', pNum, tokB{1}{1}, tokB{1}{2})));
        if ~isempty(sb)
            S = load(fullfile(behavDir, sb(1).name));
            D = S.theData; if ~isfield(S, 'theData') && isfield(S, 'data'), D = S.data; end
            if isfield(D, 'badTrial')
                [EEG, theDataGood] = dropBadTrialsAndGetGoodBehavior(EEG, D, TRIAL_START_TRIG, TRIAL_END_TRIG);
            end
        end
    end
    EEG = pop_resample(EEG, resamprate);
    EEG = pop_chanedit(EEG, 'lookup', chanLookup);
    EEG = pop_select(EEG, 'channel', 1:nChansUse);
    EEG = pop_basicfilter(EEG, 1:EEG.nbchan, 'Boundary', 'boundary', 'Cutoff', bandpass, ...
        'Design', 'butter', 'Filter', 'bandpass', 'Order', 2, 'RemoveDC', 'on');

    % Remove same channels as session 1 (to match ICA dimensionality)
    if ~isempty(badChans_s1)
        EEG = pop_select(EEG, 'nochannel', badChans_s1);
    end
    EEG = pop_reref(EEG, []);

    % Apply ICA weights from session 1
    EEG.icaweights = ICA_weights;
    EEG.icasphere  = ICA_sphere;
    EEG.icachansind = 1:EEG.nbchan;
    EEG = eeg_checkset(EEG, 'ica');

    % ICLabel & remove bad components (Muscle, Eye, Heart, Line Noise, Channel Noise, Other)
    EEG = iclabel(EEG);
    cls = EEG.etc.ic_classification.ICLabel.classifications;
    [~, maxIdx] = max(cls, [], 2);
    badIC = find(maxIdx >= 2 & maxIdx <= 6);  % 2=Muscle, 3=Eye, 4=Heart, 5=Line, 6=Channel, 7=Other
    if ~isempty(badIC)
        EEG = pop_subcomp(EEG, badIC, 0);
    end

    % Interpolate rejected channels
    EEG = pop_interp(EEG, Info.chanlocs, 'spherical');

    % Save post_ICA EEG
    outName = [basename, '_postICA'];
    pop_saveset(EEG, 'filename', [outName, '.set'], 'filepath', outDir);
    % Save good-trials-only behavior for 1:1 alignment with post_ICA
    if ~isempty(theDataGood)
        theData = theDataGood; 
        save(fullfile(behavGoodDir, [basename, '.mat']), 'theData', '-v7.3');
    end
end

%% Channel report (Oddball-style: Participant, Session, Block, RejectedChannels, RejectedList, RejectedNames)
procDir = fullfile(rootPath, 'processed', sprintf('p%d', pNum));
chanReportPath = fullfile(procDir, sprintf('p%d_ChannelReport.csv', pNum));
report = {};
report(1,:) = {'Participant','Session','Block','RejectedChannels','RejectedList','RejectedNames'};
if isempty(badChans_s1)
    rejectedList = '';
    rejectedNames = '';
else
    rejectedList = strjoin(string(badChans_s1), ',');
    rejectedNames = strjoin(labelsAll(badChans_s1), ', ');
end
for k = 1:numel(allFiles)
    tokB = regexp(allFiles(k).name, 's(\d+)b(\d+)\.set', 'tokens');
    if isempty(tokB), continue; end
    sessionNum = str2double(tokB{1}{1});
    blockNum   = str2double(tokB{1}{2});
    report(end+1,:) = {sprintf('p%d', pNum), sessionNum, blockNum, Info.nChansRejected, rejectedList, rejectedNames}; %#ok<AGROW>
end
T = cell2table(report(2:end,:), 'VariableNames', report(1,:));
writetable(T, chanReportPath);
fprintf('Channel report saved: %s\n', chanReportPath);

fprintf('Done. Post-ICA data saved to %s\n', outDir);
end

%% =====================================================================
function [EEG, theDataGood] = dropBadTrialsAndGetGoodBehavior(EEG, theData, trialStartTrig, trialEndTrig)
% Remove bad-trial segments from continuous EEG; return good-trials-only theData for alignment.
% Trial = trigger trialStartTrig (210) to trigger trialEndTrig (216). Inserts boundary events.
types = eventTypesToStr(EEG);
startIdx = find(strcmp(types, trialStartTrig));
endIdx   = find(strcmp(types, trialEndTrig));
nTrials = numel(startIdx);
if numel(endIdx) ~= nTrials
    warning('dropBadTrials: EEG has %d start (%s) vs %d end (%s) events. Skipping drop.', ...
        nTrials, trialStartTrig, numel(endIdx), trialEndTrig);
    theDataGood = theData;
    return;
end
bad = theData.badTrial(:);
if numel(bad) ~= nTrials
    warning('dropBadTrials: EEG has %d trials, behavior has %d trials. Skipping drop.', nTrials, numel(bad));
    theDataGood = theData;
    return;
end
goodIdx = find(bad == 0);
if isempty(goodIdx)
    error('dropBadTrials: no good trials in block.');
end

% Trial boundaries: trial i = startLat(i) to endLat(i)
startLat = round([EEG.event(startIdx).latency]);
endLat   = round([EEG.event(endIdx).latency]);
newData = [];
newEvents = [];
offset = 0;

for g = 1:numel(goodIdx)
    t = goodIdx(g);
    s1 = startLat(t);
    s2 = endLat(t);
    seg = EEG.data(:, s1:s2);
    newData = [newData, seg]; %#ok<AGROW>
    % Events in this segment: remap latency
    for e = 1:numel(EEG.event)
        L = EEG.event(e).latency;
        if L >= s1 && L <= s2
            ev = EEG.event(e);
            ev.latency = offset + (L - s1 + 1);
            newEvents = [newEvents, ev]; %#ok<AGROW>
        end
    end
    offset = offset + size(seg, 2);
    % Boundary between segments (except after last) - match event structure
    if g < numel(goodIdx)
        bnd = EEG.event(1);
        for fn = fieldnames(bnd)'
            bnd.(fn{1}) = [];
        end
        bnd.type = 'boundary';
        bnd.latency = offset + 0.5;
        if isfield(bnd, 'duration'), bnd.duration = 0; end
        newEvents = [newEvents, bnd]; 
    end
end

EEG.data = newData;
EEG.pnts = size(newData, 2);
EEG.xmax = EEG.xmin + (EEG.pnts - 1) / EEG.srate;
if ~isempty(newEvents)
    [~, ord] = sort([newEvents.latency]);
    EEG.event = newEvents(ord);
end
EEG = eeg_checkset(EEG, 'eventconsistency');

% Good-trials-only behavior (same order as EEG)
fn = fieldnames(theData);
theDataGood = struct();
for f = 1:numel(fn)
    v = theData.(fn{f});
    if isnumeric(v) && size(v, 1) == nTrials
        theDataGood.(fn{f}) = v(goodIdx, :);
    elseif iscell(v) && numel(v) == nTrials
        theDataGood.(fn{f}) = v(goodIdx);
    elseif isnumeric(v) && isvector(v) && numel(v) == nTrials
        theDataGood.(fn{f}) = v(goodIdx);
    else
        theDataGood.(fn{f}) = v;
    end
end
theDataGood.nTrials = numel(goodIdx);
theDataGood.badTrial = zeros(numel(goodIdx), 1); 
end

function typesStr = eventTypesToStr(EEG)
typesStr = cell(1, numel(EEG.event));
for i = 1:numel(EEG.event)
    t = EEG.event(i).type;
    if isnumeric(t), typesStr{i} = num2str(t);
    elseif isstring(t), typesStr{i} = char(t);
    else, typesStr{i} = t;
    end
end
end
