function [] = BadTrialsReport(pNum)
%% BadTrialsReport - Track bad trials per block
%
% Reads behavioral .mat files and reports bad trial counts per block.
% Uses processed/p#/matlab_files/ if available, else data/p#/matlab_files/ or p#/.

if nargin < 1
    pNum = input('Participant Number: ');
end

%% Paths
rootPath = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
procDir  = fullfile(rootPath, 'processed', sprintf('p%d', pNum), 'matlab_files');
dataDir  = fullfile(rootPath, 'data', sprintf('p%d', pNum), 'matlab_files');
fallback = fullfile(rootPath, sprintf('p%d', pNum));

if exist(procDir, 'dir')
    behavDir = procDir;
    src = 'processed';
elseif exist(dataDir, 'dir')
    behavDir = dataDir;
    src = 'data';
elseif exist(fallback, 'dir')
    behavDir = fallback;
    src = 'participant folder';
else
    error('No behavior folder found for participant %d', pNum);
end

files = dir(fullfile(behavDir, sprintf('p%ds*b*.mat', pNum)));
if isempty(files)
    error('No behavioral .mat files in %s', behavDir);
end

% Sort by block number
tok = regexp({files.name}, 'p\d+s\d+b(\d+)\.mat', 'tokens');
bnum = nan(numel(files), 1);
for k = 1:numel(files)
    if ~isempty(tok{k})
        bnum(k) = str2double(tok{k}{1}{1});
    end
end
[~, ord] = sort(bnum);
files = files(ord);
bnum = bnum(ord);

%% Collect stats per block
nBlocks = numel(files);
blockNum = nan(nBlocks, 1);
nTrials  = nan(nBlocks, 1);
nBad     = nan(nBlocks, 1);
pctBad   = nan(nBlocks, 1);

for k = 1:nBlocks
    S = load(fullfile(behavDir, files(k).name));
    if isfield(S, 'theData'), D = S.theData; elseif isfield(S, 'data'), D = S.data; else, continue; end
    if ~isfield(D, 'badTrial')
        continue;
    end
    bad = D.badTrial(:);
    n = numel(bad);
    nBad(k) = sum(bad > 0);
    nTrials(k) = n;
    pctBad(k) = 100 * nBad(k) / max(1, n);
    blockNum(k) = bnum(k);
end

%% Print table
fprintf('\n=== Bad Trials Report: Participant %d (source: %s) ===\n\n', pNum, src);
fprintf('Block\tTrials\tBad\t%% Bad\n');
fprintf('-----\t------\t---\t-----\n');
for k = 1:nBlocks
    fprintf('%d\t%d\t%d\t%.1f\n', blockNum(k), nTrials(k), nBad(k), pctBad(k));
end
fprintf('-----\t------\t---\t-----\n');
fprintf('Total\t%d\t%d\t%.1f\n\n', sum(nTrials, 'omitnan'), sum(nBad, 'omitnan'), ...
    100 * sum(nBad, 'omitnan') / max(1, sum(nTrials, 'omitnan')));

%% Save Excel
outFile = fullfile(rootPath, 'processed', sprintf('p%d', pNum), sprintf('p%d_BadTrialsReport.xlsx', pNum));
outDir = fileparts(outFile);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
T = table(blockNum, nTrials, nBad, pctBad, 'VariableNames', {'Block', 'nTrials', 'nBad', 'pctBad'});
writetable(T, outFile);
fprintf('Saved: %s\n', outFile);
end
