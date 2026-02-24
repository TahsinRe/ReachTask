function [] = ChannelReport(pNum)
%% ChannelReport - Add channel report sheet to BadTrialsReport.xlsx
%
% Reads p{N}_ChannelReport.csv from Preprocess_ReachTask (Oddball-style:
% Participant, Session, Block, RejectedChannels, RejectedList, RejectedNames),
% and adds it as sheet "ChannelReport" to BadTrialsReport.
%
% Run after BadTrialsReport and Preprocess_ReachTask.

if nargin < 1
    pNum = input('Participant Number: ');
end

rootPath = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
procDir  = fullfile(rootPath, 'processed', sprintf('p%d', pNum));
outFile  = fullfile(procDir, sprintf('p%d_BadTrialsReport.xlsx', pNum));
chanReportPath = fullfile(procDir, sprintf('p%d_ChannelReport.csv', pNum));

if ~exist(outFile, 'file')
    error('BadTrialsReport Excel not found: %s (run BadTrialsReport first)', outFile);
end

if ~exist(chanReportPath, 'file')
    error('Channel report not found: %s\nRun Preprocess_ReachTask first.', chanReportPath);
end

Tchan = readtable(chanReportPath);
writetable(Tchan, outFile, 'Sheet', 'ChannelReport');

fprintf('Added ChannelReport sheet to %s\n', outFile);

end
