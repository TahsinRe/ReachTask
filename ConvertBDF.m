function [] = ConvertBDF(pNum)
% Converts .bdf files (block-by-block) into EEGLAB .set/.fdt files.
if nargin < 1
    pNum = input('Participant Number: ');
end
%% --- Open EEGLAB ---
cd('C:/Users/tahsi/OneDrive/Desktop/EEG/eeglab2024.0'); eeglab; close;

%% --- Paths ---
rootPath = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
rawPath  = [rootPath, '/data/p', num2str(pNum), '/raw/'];
savePath = [rootPath, '/processed/p', num2str(pNum), '/eeglab/'];
nBlocks     = 32;
maxSessions = 1;   % sessions 1–3 

%% --- Session and block loop ---
for sessionNum = 1:maxSessions

    for iBlock = 1:nBlocks

        % Raw .bdf file
        bdfFile = [rawPath, 'p', num2str(pNum), 's', num2str(sessionNum), 'b', num2str(iBlock), '.bdf'];
        % Output name
        outName = ['p', num2str(pNum), 's', num2str(sessionNum), 'b', num2str(iBlock)];

        % If .bdf doesn't exist → silently skip
        if exist(bdfFile, 'file') ~= 2
            continue;
        end

        % If .set exists already → silently skip
        if exist([savePath, outName, '.set'], 'file') == 2
            continue;
        end

        % Create output directory if missing
        if ~exist(savePath, 'dir')
            mkdir(savePath);
        end

        % Import BDF and save SET/FDT
        EEG = pop_biosig(bdfFile);

        pop_saveset(EEG, 'filename', [outName, '.set'], 'filepath', savePath);

    end
end

% Comment out the next line if you only want conversion and no preprocessing:
%Stage2_Preprocess(pNum);

end
