function [] = Beta_Weights_ReachTask_AzEl(pNum, epochTypes, theoreticalPath, doFig, swapAzEl)
%% Beta_Weights_ReachTask_AzEl - Regression of data RDM on 3 model RDMs (no GoColor)
%
% Three models: (1) Azimuth, (2) Elevation, (3) Location (az+el).
% Location = same vs different (0/1). Note: in 24-cond design Location has no variance (all pairs differ).
% vec_neural ~ [1, Az, El, Location]. Beta 2-4 = Az, El, Location.
%
% Loads: results/RSA_AzEl/RDA_p<N>_<Epoch>.mat (24 cond), theoretical_RDM_ReachTask_AzEl.mat.
% Saves: results/RSA_AzEl/beta_weights_p<N>_<Epoch>.mat

if nargin < 1
    pNum = input('Participant Number: ');
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
if nargin < 5 || isempty(swapAzEl)
    swapAzEl = false;  % set true if Az/El labels appear swapped (e.g. Az flat, El has dynamics)
end

if ~exist(theoreticalPath, 'file')
    error('Theoretical RDMs not found: %s. Run Theoretical_RDM_ReachTask_AzEl first.', theoreticalPath);
end
load(theoreticalPath, 'model_azimuth', 'model_elevation', 'model_integration', 'nC', 'azimuth_by_loc', 'elevation_by_loc');

rootPath = 'C:/Users/tahsi/OneDrive/Desktop/EEG/Experiments/ReachTask';
rdaDir = fullfile(rootPath, 'results', 'RSA_AzEl');

idx = find(tril(ones(nC), 0));  % include diagonal (same vs same = 0)
if swapAzEl
    vec_azimuth   = model_elevation(idx);  % swap: display "Az" = elevation
    vec_elevation = model_azimuth(idx);   % display "El" = azimuth
else
    vec_azimuth   = model_azimuth(idx);
    vec_elevation = model_elevation(idx);
end
vec_integration = model_integration(idx);
X = [ones(length(idx), 1), vec_azimuth, vec_elevation, vec_integration];

nPred = 3;  % Az, El, Location
beta_by_epoch = cell(numel(epochTypes), 1);
times_by_epoch = cell(numel(epochTypes), 1);

for e = 1:numel(epochTypes)
    epName = epochTypes{e};
    rdaFile = fullfile(rdaDir, sprintf('RDA_p%d_%s.mat', pNum, epName));
    if ~exist(rdaFile, 'file')
        warning('Skip %s: no RDA file %s', epName, rdaFile);
        continue;
    end
    load(rdaFile, 'RDM_all', 'times_bins', 'nTimeBins');
    if isempty(times_bins)
        times_bins = 1:nTimeBins;
    end

    beta_all = zeros(nPred + 1, nTimeBins);

    for t = 1:nTimeBins
        neural_RDM = RDM_all(:, :, t);
        neural_RDM(1:nC+1:end) = 0;  % diagonal = 0 (same vs same; RDA leaves it NaN)
        vec_neural = neural_RDM(idx);
        valid = isfinite(vec_neural);
        if sum(valid) < nPred + 5
            beta_all(:, t) = NaN;
            continue;
        end
        y = vec_neural(valid);
        Xv = X(valid, :);
        beta_all(:, t) = regress(y, Xv);
    end

    outFile = fullfile(rdaDir, sprintf('beta_weights_p%d_%s.mat', pNum, epName));
    save(outFile, 'beta_all', 'times_bins', 'nTimeBins', 'nPred', 'X');
    fprintf('Saved %s (beta 2-4: Az, El, Location)\n', outFile);

    beta_by_epoch{e} = beta_all;
    times_by_epoch{e} = times_bins;

    if doFig
        predNames = {'Azimuth', 'Elevation', 'Integrated'};
        epToXLab = struct('Target','Time Relative to Preview (ms)', 'GoCue','Time Relative to Go Cue (ms)', ...
            'Movement','Time Relative to Movement (ms)', 'Touch','Time Relative to Touch (ms)');
        if isfield(epToXLab, epName), xLab = epToXLab.(epName); else, xLab = sprintf('Time relative to %s (ms)', epName); end
        colors = [0.2 0.4 0.8; 0.9 0.4 0.6; 0.2 0.6 0.35];  % Az=blue, El=pink, Loc=green
        tms = times_bins;
        fig = figure('Visible', 'off', 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial');
        for p = 1:3
            b = beta_all(p+1, :);
            valid = isfinite(b);
            if any(valid)
                yRow = b(valid);
                pad = max(0.1, (max(yRow)-min(yRow))*0.05);
                yLims = [min(yRow)-pad, max(yRow)+pad];
            else
                yLims = [-0.5 0.5];
            end
            ax = subplot(3, 1, p);
            b = beta_all(p+1, :);
            valid = isfinite(b);
            if any(valid)
                plot(ax, tms, b, 'Color', colors(p,:), 'LineWidth', 1.2);
            end
            hold(ax, 'on');
            plot(ax, [min(tms), max(tms)], [0 0], 'k-', 'LineWidth', 0.5);
            plot(ax, [0 0], yLims, 'k--', 'LineWidth', 0.5);
            ylabel(ax, predNames{p}, 'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 9);
            set(ax, 'XLim', [min(tms), max(tms)], 'YLim', yLims, 'Box', 'off', 'FontName', 'Arial');
            grid(ax, 'off');
            if p == 3, xlabel(ax, xLab, 'FontName', 'Arial', 'FontSize', 9); end
        end
        figPath = fullfile(rdaDir, sprintf('beta_weights_p%d_%s.png', pNum, epName));
        saveas(fig, figPath);
        saveas(fig, fullfile(rdaDir, sprintf('beta_weights_p%d_%s.fig', pNum, epName)));
        close(fig);
        fprintf('  Saved figure: %s\n', figPath);
    end
end

% Combined figure: theoretical RDMs in left column, 4 epochs as columns, 3 predictors as rows (5x3 layout)
if doFig && numel(epochTypes) >= 1 && any(~cellfun(@isempty, beta_by_epoch))
    fig = figure('Visible', 'off', 'Position', [50 50 1400 700], 'Color', 'w', 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial');
    predNames = {'Azimuth', 'Elevation', 'Integrated'};
    xLabels = {'Time Relative to Preview (ms)', 'Time Relative to Go Cue (ms)', ...
               'Time Relative to Movement (ms)', 'Time Relative to Touch (ms)'};
    colors = [0.2 0.4 0.8; 0.9 0.4 0.6; 0.2 0.6 0.35];  % Az=blue, El=pink, Loc=green
    nE = numel(epochTypes);
    gapLabelPlot = 0.025;   % distance from label to plot edge (gapH must be >= gapLabelPlot)
    marginL = 0.08; marginR = 0.02; marginB = 0.08 + 2*gapLabelPlot; marginT = 0.06; gapH = 0.04; gapV = 0.04;  % marginB + 2*gapLabelPlot: x-axis label same distance as Regression Beta Weights
    wRDM = 0.16;   % width of theoretical RDM column (bigger model plots)
    w = (1 - marginL - marginR - wRDM - gapH - 3*gapH) / 4;   % 4 columns (epochs), after RDM column
    h = (1 - marginB - marginT - 2*gapV) / 3;   % 3 rows (predictors)
    % Reorder RDMs by (azimuth, elevation) for display
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
    % First pass: get y-axis limits per row (predictor) so each row scales to show variation
    yLimsByRow = cell(3, 1);
    for p = 1:3
        yRow = [];
        for ee = 1:nE
            if isempty(beta_by_epoch{ee}), continue; end
            bet = beta_by_epoch{ee};
            b = bet(p+1, :);
            valid = isfinite(b);
            if any(valid), yRow = [yRow, b(valid)]; end
        end
        if isempty(yRow)
            yLimsByRow{p} = [-0.5 0.5];
        else
            pad = max(0.1, (max(yRow)-min(yRow))*0.05);
            yLimsByRow{p} = [min(yRow)-pad, max(yRow)+pad];
        end
    end
    % First: plot theoretical RDMs in left column (all same size, no colorbar)
    xRDM = marginL;
    for p = 1:3
        axRDM = axes('Position', [xRDM, 1 - marginT - p*(h+gapV), wRDM*0.9, h]);
        imagesc(axRDM, rdmMats{p});
        set(axRDM, 'Color', 'w', 'YDir', 'normal', 'XLim', [0.5 nC+0.5], 'YLim', [0.5 nC+0.5], 'Box', 'off', 'FontName', 'Arial', 'FontSize', 7);
        colormap(axRDM, cmap);
        caxis(axRDM, [0 1]);
        axis(axRDM, 'image');
        if p == 1, title(axRDM, 'Model', 'FontSize', 9, 'FontWeight', 'bold', 'FontName', 'Arial'); end
    end
    % Second: plot beta weights (4 event columns to the right of RDMs)
    xPlotStart = marginL + wRDM + gapH;
    for p = 1:3
        yLims = yLimsByRow{p};
        for ee = 1:nE
            if isempty(beta_by_epoch{ee}), continue; end
            epName = epochTypes{ee};
            tms = times_by_epoch{ee};
            bet = beta_by_epoch{ee};
            col = ee; row = p;
            ax = axes('Position', [xPlotStart + (col-1)*(w+gapH), 1 - marginT - row*(h+gapV), w, h]);
            b = bet(p+1, :);
            valid = isfinite(b);
            if any(valid)
                plot(ax, tms, b, 'Color', colors(p,:), 'LineWidth', 1.3);
            end
            hold(ax, 'on');
            plot(ax, [min(tms), max(tms)], [0 0], 'k-', 'LineWidth', 0.5);
            plot(ax, [0 0], yLims, 'k--', 'LineWidth', 0.5);
            if p == 1, title(ax, epName, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial'); end
            if p == 3
                xlh = xlabel(ax, xLabels{ee}, 'FontSize', 9, 'FontName', 'Arial');
                xlh.Units = 'normalized';
                xlh.Position(2) = xlh.Position(2) - 0.04;   % more gap between tick numbers and "Time Relative to..." text
            end
            set(ax, 'XLim', [min(tms), max(tms)], 'YLim', yLims, 'Box', 'off', 'FontSize', 8, 'FontName', 'Arial');
            grid(ax, 'off');
        end
    end
    % Add labels: Azimuth, Elevation, Integrated first, then RDMs, then Regression Beta Weights (after theoretical model)
    axLabel = axes('Position', [0 0 1 1], 'Visible', 'off');
    xPred = 0.03;       % Azimuth, Elevation, Integrated (closer to Model column)
    xYLabel = marginL + wRDM + gapH - gapLabelPlot;  % Regression Beta Weights: gapLabelPlot from plot edge
    for p = 1:3
        yRow = 1 - marginT - (p-0.5)*(h+gapV);
        text(axLabel, xPred, yRow, predNames{p}, 'FontSize', 9, 'FontWeight', 'bold', 'FontName', 'Arial', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90, 'Units', 'normalized');
    end
    text(axLabel, xYLabel, 0.5, 'Regression Beta Weights', 'FontSize', 9, 'FontName', 'Arial', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90, 'Units', 'normalized');
    figPath = fullfile(rdaDir, sprintf('beta_weights_p%d_allEpochs_stacked.png', pNum));
    print(fig, figPath, '-dpng', '-r300');  % 300 dpi for publication
    saveas(fig, fullfile(rdaDir, sprintf('beta_weights_p%d_allEpochs_stacked.fig', pNum)));
    close(fig);
    fprintf('  Saved combined figure: %s\n', figPath);
end

fprintf('\nBeta weights (AzEl only) done. 3 models: beta_all(2)=Az, (3)=El, (4)=Location.\n');
end
