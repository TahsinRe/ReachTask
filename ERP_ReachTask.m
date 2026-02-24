function [] = ERP_ReachTask(pNum)
%% ERP_ReachTask - Create epochs from post_ICA data (like ERP_Kirtan)
%
% Loads post_ICA .set files and creates four epoch types:
%   1. Target (211):  -100 to 1000 ms
%   2. GoCue (213):   -100 to 100 ms
%   3. Movement (214): -100 to 500 ms
%   4. Touch (215):   -100 to 300 ms (shorter to avoid truncation)
%
% Baseline correction: subtract mean of -100 to 0 ms (pre-event) per channel per epoch.
% Run after Preprocess_ReachTask. Saves to ERPs_Target/, ERPs_GoCue/, etc.

if nargin < 1
    pNum = input('Participant Number: ');
end

%% --- Open EEGLAB ---
cd('C:/Users/tahsi/OneDrive/Desktop/EEG/eeglab2024.0'); eeglab; close;

%% Paths
rootPath = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
inDir   = fullfile(rootPath, 'processed', sprintf('p%d', pNum), 'post_ICA');

if ~exist(inDir, 'dir')
    error('post_ICA not found: %s (run Preprocess_ReachTask first)', inDir);
end

%% Epoch definitions: {name, trigger, timelim_sec [pre, post]}
epochDefs = {
    'Target',   '211', [-0.1  1.0];   % -100 to 1000 ms
    'GoCue',    '213', [-0.1  0.1];   % -100 to 100 ms
    'Movement', '214', [-0.1  0.5];   % -100 to 500 ms
    'Touch',    '215', [-0.1  0.3]    % -100 to 300 ms 
};

%% Get post_ICA files (sorted by session, block)
files = dir(fullfile(inDir, 'p*_postICA.set'));
if isempty(files)
    error('No post_ICA .set files in %s', inDir);
end
tok = regexp({files.name}, 'p\d+s(\d+)b(\d+)_postICA', 'tokens');
valid = ~cellfun(@isempty, tok);
files = files(valid);
snums = cellfun(@(x) str2double(x{1}{1}), tok(valid));
bnums = cellfun(@(x) str2double(x{1}{2}), tok(valid));
[~, ord] = sortrows([snums(:), bnums(:)], [1 2]);
files = files(ord);

%% Loop over epoch types 
for eType = 1:size(epochDefs, 1)
    epName   = epochDefs{eType, 1};
    trigger  = epochDefs{eType, 2};
    timelim  = epochDefs{eType, 3};
    savepath = fullfile(rootPath, 'processed', sprintf('p%d', pNum), ['ERPs_', epName]);

    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end

    fprintf('\n=== Epoch type: %s (trigger %s, [%.1f %.1f] s) ===\n', epName, trigger, timelim(1), timelim(2));

    for fIdx = 1:numel(files)
        fname = files(fIdx).name;
        fprintf('  %s\n', fname);

        EEG = pop_loadset(fname, inDir);

        % Epoch around trigger
        EEG = pop_epoch(EEG, {trigger}, timelim, 'newname', [epName '_epoch'], 'epochinfo', 'yes');

        % Baseline correction: -100 to 0 ms (pre-event) per channel per epoch
        baseline_ms = [timelim(1)*1000, 0];
        EEG = pop_rmbase(EEG, baseline_ms);

        % Reject epochs with extreme values (optional; remove if not desired)
        % EEG = pop_autorej(EEG, 'nogui', 'on');

        % Save epoched data
        [~, base, ~] = fileparts(fname);
        base = strrep(base, '_postICA', '');
        outName = [base, '_', epName, 'Epoch.set'];
        pop_saveset(EEG, 'filename', outName, 'filepath', savepath);
    end
end

fprintf('\nDone. Epochs saved to ERPs_Target, ERPs_GoCue, ERPs_Movement, ERPs_Touch.\n');


end
