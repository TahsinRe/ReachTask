function [] = Beta_Weights_Group_ReachTask_AzEl(pList, epochTypes, theoreticalPath, doFig)
%% Beta_Weights_Group_ReachTask_AzEl - Group mean of regression beta weights across participants
%
% Loads beta_weights_p<N>_<Epoch>.mat for each participant, computes mean and SEM,
% plots with shaded p-variation (SEM) region. Same layout as individual Beta_Weights_ReachTask_AzEl.
%
% Inputs:
%   pList       - vector of participant numbers (e.g. [1 5 6])
%   epochTypes  - cell of epoch names (default: Target, GoCue, Movement, Touch)
%   theoreticalPath - path to theoretical_RDM_ReachTask_AzEl.mat
%   doFig       - if true, create and save figure (default: true)
%
% Saves: results/RSA_AzEl/beta_weights_group_allEpochs_stacked.png
%        results/RSA_AzEl/beta_weights_group_allEpochs_stacked.mat (mean, SEM, times)

if nargin < 1 || isempty(pList)
    pList = input('Participant list (e.g. [1 5 6]): ');
end
if nargin < 2 || isempty(epochTypes)
    epochTypes = {'Target','GoCue','Movement','Touch'};
end
if nargin < 3 || isempty(theoreticalPath)
    rootPath = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
    theoreticalPath = fullfile(rootPath, 'results', 'RSA_AzEl', 'theoretical_RDM_ReachTask_AzEl.mat');
end
if nargin < 4 || isempty(doFig)
    doFig = true;
end

if ~exist(theoreticalPath, 'file')
    error('Theoretical RDMs not found: %s. Run Theoretical_RDM_ReachTask_AzEl first.', theoreticalPath);
end
load(theoreticalPath, 'model_azimuth', 'model_elevation', 'model_integration', 'nC', 'azimuth_by_loc', 'elevation_by_loc');

rootPath = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
rdaDir = fullfile(rootPath, 'results', 'RSA_AzEl');

nE = numel(epochTypes);
beta_mean_by_epoch = cell(nE, 1);
beta_sem_by_epoch  = cell(nE, 1);
times_by_epoch      = cell(nE, 1);

for e = 1:nE
    epName = epochTypes{e};
    betaCell = {};
    tmsRef = [];
    for ip = 1:numel(pList)
        pNum = pList(ip);
        bwFile = fullfile(rdaDir, sprintf('beta_weights_p%d_%s.mat', pNum, epName));
        if ~exist(bwFile, 'file')
            warning('Skip p%d %s: no file %s', pNum, epName, bwFile);
            continue;
        end
        load(bwFile, 'beta_all', 'times_bins');
        if isempty(tmsRef), tmsRef = times_bins; end
        betaCell{end+1} = beta_all(2:4, :)';  % [nTimeBins x 3]
    end
    if isempty(betaCell)
        continue;
    end
    % Align times: use shortest nTimeBins across participants
    nT = size(betaCell{1}, 1);
    for i = 2:numel(betaCell)
        nT = min(nT, size(betaCell{i}, 1));
    end
    B = nan(numel(betaCell), nT, 3);
    for i = 1:numel(betaCell)
        B(i, :, :) = betaCell{i}(1:nT, :);
    end
    beta_mean = squeeze(nanmean(B, 1));  % [nT x 3]
    beta_sem  = squeeze(nanstd(B, 0, 1) / sqrt(size(B, 1)));  % [nT x 3]
    beta_mean_by_epoch{e} = beta_mean';
    beta_sem_by_epoch{e}  = beta_sem';
    times_by_epoch{e} = tmsRef(1:nT);
end

% Save group data
outMat = fullfile(rdaDir, 'beta_weights_group_allEpochs_stacked.mat');
save(outMat, 'beta_mean_by_epoch', 'beta_sem_by_epoch', 'times_by_epoch', 'pList', 'epochTypes');
fprintf('Saved group beta weights: %s (n=%d participants)\n', outMat, numel(pList));

% Figure: same layout as individual, with mean + shaded SEM
if ~doFig || all(cellfun(@isempty, beta_mean_by_epoch))
    return;
end

predNames = {'Azimuth', 'Elevation', 'Integrated'};
xLabels = {'Time Relative to Preview (ms)', 'Time Relative to Go Cue (ms)', ...
           'Time Relative to Movement (ms)', 'Time Relative to Touch (ms)'};
colors = [0.6 0.2 0.9; 1 0.3 0.6; 0.2 0.75 0.5];  % Az=purple, El=pink, Int=teal-green
gapLabelPlot = 0.025;
marginL = 0.08; marginR = 0.02; marginB = 0.08 + 2*gapLabelPlot; marginT = 0.06; gapH = 0.04; gapV = 0.04;
wRDM = 0.16;
w = (1 - marginL - marginR - wRDM - gapH - 3*gapH) / 4;
h = (1 - marginB - marginT - 2*gapV) / 3;

az_vec = azimuth_by_loc(:);
el_vec = elevation_by_loc(:);
[~, ord] = sortrows([az_vec, el_vec], [1 2]);
M_az = model_azimuth(ord, ord);
M_el = model_elevation(ord, ord);
M_int = model_integration(ord, ord);
rdmMats = {M_az, M_el, M_int};
blue = [50, 102, 152] / 255;
red  = [152, 0, 49] / 255;
cmap = interp1([0 0.5 1], [blue; 1 1 1; red], linspace(0, 1, 256));

% Y limits from mean ± SEM
yLimsByRow = cell(3, 1);
for p = 1:3
    yRow = [];
    for ee = 1:nE
        if isempty(beta_mean_by_epoch{ee}), continue; end
        m = beta_mean_by_epoch{ee};
        s = beta_sem_by_epoch{ee};
        b = m(p, :);
        v = isfinite(b);
        if any(v)
            yRow = [yRow, b(v), (b(v) + s(p,v)), (b(v) - s(p,v))];
        end
    end
    if isempty(yRow)
        yLimsByRow{p} = [-0.5 0.5];
    else
        pad = max(0.1, (max(yRow)-min(yRow))*0.05);
        yLimsByRow{p} = [min(yRow)-pad, max(yRow)+pad];
    end
end

fig = figure('Visible', 'off', 'Position', [50 50 1400 700], 'Color', 'w', 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial');

% Theoretical RDMs
xRDM = marginL;
for p = 1:3
    axRDM = axes('Position', [xRDM, 1 - marginT - p*(h+gapV), wRDM*0.9, h]);
    imagesc(axRDM, rdmMats{p});
    set(axRDM, 'Color', 'w', 'YDir', 'normal', 'XLim', [0.5 nC+0.5], 'YLim', [0.5 nC+0.5], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 7);
    colormap(axRDM, cmap);
    caxis(axRDM, [0 1]);
    axis(axRDM, 'image');
end

% Beta weight plots with shaded SEM
xPlotStart = marginL + wRDM + gapH;
for p = 1:3
    yLims = yLimsByRow{p};
    for ee = 1:nE
        if isempty(beta_mean_by_epoch{ee}), continue; end
        epName = epochTypes{ee};
        tms = times_by_epoch{ee};
        m = beta_mean_by_epoch{ee};
        s = beta_sem_by_epoch{ee};
        col = ee; row = p;
        ax = axes('Position', [xPlotStart + (col-1)*(w+gapH), 1 - marginT - row*(h+gapV), w, h]);
        b = m(p, :);
        sem = s(p, :);
        valid = isfinite(b);
        if any(valid)
            % Shaded region: mean ± SEM (p-variation)
            fillX = [tms(valid), fliplr(tms(valid))];
            fillY = [b(valid) + sem(valid), fliplr(b(valid) - sem(valid))];
            fill(ax, fillX, fillY, colors(p,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            hold(ax, 'on');
            plot(ax, tms, b, 'Color', colors(p,:), 'LineWidth', 1.3);
        end
        hold(ax, 'on');
        plot(ax, [min(tms), max(tms)], [0 0], 'k-', 'LineWidth', 0.5);
        plot(ax, [0 0], yLims, 'k--', 'LineWidth', 0.5);
        if p == 1, title(ax, epName, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial'); end
        if p == 3
            xlh = xlabel(ax, xLabels{ee}, 'FontSize', 9, 'FontName', 'Arial');
            xlh.Units = 'normalized';
            xlh.Position(2) = xlh.Position(2) - 0.04;
        end
        set(ax, 'XLim', [min(tms), max(tms)], 'YLim', yLims, 'Box', 'off', 'FontSize', 8, 'FontName', 'Arial');
        grid(ax, 'off');
    end
end

% Labels: Azimuth, Elevation, Integrated + Regression Beta Weights
axLabel = axes('Position', [0 0 1 1], 'Visible', 'off');
xPred = 0.06;   % closer to Model RDMs
xYLabel = marginL + wRDM + gapH - gapLabelPlot;
for p = 1:3
    yRow = 1 - marginT - (p-0.5)*(h+gapV);
    text(axLabel, xPred, yRow, predNames{p}, 'FontSize', 9, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90, 'Units', 'normalized');
end
text(axLabel, xYLabel, 0.5, 'Regression Beta Weights', 'FontSize', 9, 'FontName', 'Arial', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90, 'Units', 'normalized');

figPath = fullfile(rdaDir, 'beta_weights_group_allEpochs_stacked.png');
print(fig, figPath, '-dpng', '-r300');
saveas(fig, fullfile(rdaDir, 'beta_weights_group_allEpochs_stacked.fig'));
close(fig);
fprintf('  Saved group figure: %s\n', figPath);

fprintf('\nGroup beta weights done. n=%d participants.\n', numel(pList));
