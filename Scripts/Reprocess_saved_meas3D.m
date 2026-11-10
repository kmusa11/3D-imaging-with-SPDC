%% Reconstruct 3D depth images from measured SPDC-SFG delay scans

clear; clc; close all;

[fileName, filePath] = uigetfile('*.mat', 'Select a saved 3D measurement');
if isequal(fileName, 0)
    return;
end

S = load(fullfile(filePath, fileName));
meas3D = double(S.meas3D);
xs = S.xs(:).';
ys = S.ys(:).';
delayRange = S.delayRange(:).';

% Reconstruction settings
h = 0.995;
thrCV = 0.015;
minMean = 0;

showDepthMap = true;
showClusterPlot = false;
saveProcessedData = false;

fontSize_axes = 20;
fontSize_labels = 24;
zLimits_um = [-150, 520];

[Ny, Nx, nDelay] = size(meas3D);
if numel(xs) ~= Nx || numel(ys) ~= Ny || numel(delayRange) ~= nDelay
    error('The meas3D dimensions do not match xs, ys, and delayRange.');
end
if any(diff(delayRange) <= 0)
    error('delayRange must be strictly increasing.');
end
if h <= 0 || h >= 1 || thrCV < 0
    error('Use a crossing level between 0 and 1 and a nonnegative CV threshold.');
end

centerDelayMap_ps = nan(Ny, Nx);
leftCrossing_ps = nan(Ny, Nx);
rightCrossing_ps = nan(Ny, Nx);
cvMap = nan(Ny, Nx);
meanMap = nan(Ny, Nx);
maxVals = nan(Ny, Nx);

for iy = 1:Ny
    for ix = 1:Nx
        trace = squeeze(meas3D(iy, ix, :)).';
        if any(~isfinite(trace)) || max(trace) <= 0
            continue;
        end

        m = mean(trace);
        if m <= 0
            continue;
        end

        meanMap(iy, ix) = m;
        cvMap(iy, ix) = std(trace, 0) / abs(m);
        maxVals(iy, ix) = max(trace);

        [centerDelayMap_ps(iy, ix), leftCrossing_ps(iy, ix), ...
            rightCrossing_ps(iy, ix)] = crossingCenter(delayRange, trace, h);
    end
end

validPix = isfinite(centerDelayMap_ps) & isfinite(cvMap) & ...
    cvMap > thrCV & meanMap > minMean;
if ~any(validPix(:))
    error('No valid pixels for h = %.4f and CV threshold = %.4f.', h, thrCV);
end

refDelay_ps = median(centerDelayMap_ps(validPix));
centerDepthMap_um = (centerDelayMap_ps - refDelay_ps) * (299.792458 / 2);
centerDepthMap_um(~validPix) = NaN;

maxValsCV = maxVals;
maxValsCV(~validPix) = NaN;
normalizedCounts = maxValsCV / max(maxValsCV(validPix));

fprintf('File: %s\n', fileName);
fprintf('Data size: %d x %d x %d\n', Ny, Nx, nDelay);
fprintf('h = %.4f, CV threshold = %.4f\n', h, thrCV);
fprintf('Accepted pixels: %d / %d\n', nnz(validPix), Nx * Ny);

[X2, Y2] = meshgrid(xs, ys);

figure('Color', 'w', 'Name', '3D depth reconstruction', ...
    'Position', [200 120 900 720]);
ax = axes;
scatter3(ax, X2(validPix), Y2(validPix), centerDepthMap_um(validPix), ...
    40, normalizedCounts(validPix), 'filled');
xlabel(ax, '$x$ position', 'Interpreter', 'latex', 'FontSize', fontSize_labels);
ylabel(ax, '$y$ position', 'Interpreter', 'latex', 'FontSize', fontSize_labels);
zlabel(ax, 'Relative depth ($\mu$m)', 'Interpreter', 'latex', 'FontSize', fontSize_labels);
xlim(ax, [min(xs), max(xs)]);
ylim(ax, [min(ys), max(ys)]);
zlim(ax, zLimits_um);
clim(ax, [0, 1]);
colorbar(ax);
ax.Position = [0.2 0.17 0.58 0.72];
view(ax, 3);
pbaspect(ax, [1 1 0.75]);
axis(ax, 'vis3d');
grid(ax, 'on');
box(ax, 'on');
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fontSize_axes, ...
    'LineWidth', 1.2, 'TickDir', 'out');

if showDepthMap
    figure('Color', 'w', 'Name', 'Relative depth map');
    imagesc(xs, ys, centerDepthMap_um);
    axis xy;
    axis image;
    xlabel('$x$ position', 'Interpreter', 'latex', 'FontSize', fontSize_labels);
    ylabel('$y$ position', 'Interpreter', 'latex', 'FontSize', fontSize_labels);
    cb = colorbar;
    ylabel(cb, 'Relative depth ($\mu$m)', 'Interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontSize_axes);
end

depthValues = centerDepthMap_um(validPix);
if numel(depthValues) >= 4
    [clusterID, clusterMeans] = kmeans1D_twoClusters(depthValues);
    [~, order] = sort(clusterMeans);
    lowerSurface = depthValues(clusterID == order(1));
    upperSurface = depthValues(clusterID == order(2));

    lowMean_um = mean(lowerSurface);
    highMean_um = mean(upperSurface);
    separation_um = highMean_um - lowMean_um;
    separationSEM_um = sqrt(std(lowerSurface)^2 / numel(lowerSurface) + ...
        std(upperSurface)^2 / numel(upperSurface));

    fprintf('Surface separation: %.2f +/- %.2f um\n', ...
        separation_um, separationSEM_um);

    if showClusterPlot
        clusterMap = nan(Ny, Nx);
        clusterMap(validPix) = clusterID;
        figure('Color', 'w', 'Name', 'Two-surface reconstruction');
        hold on;
        maskLow = clusterMap == order(1);
        maskHigh = clusterMap == order(2);
        scatter3(X2(maskLow), Y2(maskLow), centerDepthMap_um(maskLow), ...
            40, 'b', 'filled');
        scatter3(X2(maskHigh), Y2(maskHigh), centerDepthMap_um(maskHigh), ...
            40, 'r', 'filled');
        xlabel('$x$ position', 'Interpreter', 'latex');
        ylabel('$y$ position', 'Interpreter', 'latex');
        zlabel('Relative depth ($\mu$m)', 'Interpreter', 'latex');
        grid on; box on; view(3); pbaspect([1 1 0.75]);
    end
else
    separation_um = NaN;
    separationSEM_um = NaN;
end

if saveProcessedData
    [~, baseName] = fileparts(fileName);
    outputName = fullfile(filePath, ['processed_' baseName '.mat']);
    save(outputName, 'centerDelayMap_ps', 'centerDepthMap_um', ...
        'leftCrossing_ps', 'rightCrossing_ps', 'cvMap', 'meanMap', ...
        'maxVals', 'validPix', 'refDelay_ps', 'separation_um', ...
        'separationSEM_um', 'h', 'thrCV', 'minMean');
end

function [center, left, right] = crossingCenter(t, trace, h)
    center = NaN; left = NaN; right = NaN;
    trace = trace / max(trace);
    [~, peak] = max(trace);

    leftIdx = find(trace(1:peak) <= h, 1, 'last');
    rightOffset = find(trace(peak:end) <= h, 1, 'first');
    if isempty(leftIdx) || leftIdx >= peak || ...
            isempty(rightOffset) || rightOffset <= 1
        return;
    end

    rightIdx = peak + rightOffset - 1;
    if trace(leftIdx + 1) == trace(leftIdx) || ...
            trace(rightIdx) == trace(rightIdx - 1)
        return;
    end

    left = t(leftIdx) + (h - trace(leftIdx)) * ...
        (t(leftIdx + 1) - t(leftIdx)) / (trace(leftIdx + 1) - trace(leftIdx));
    right = t(rightIdx - 1) + (h - trace(rightIdx - 1)) * ...
        (t(rightIdx) - t(rightIdx - 1)) / (trace(rightIdx) - trace(rightIdx - 1));
    center = (left + right) / 2;
end

function [groups, means] = kmeans1D_twoClusters(x)
    x = x(:);
    sortedX = sort(x);
    c1 = interpPercentile(sortedX, 25);
    c2 = interpPercentile(sortedX, 75);
    if c1 == c2
        c1 = min(x);
        c2 = max(x);
    end
    if c1 == c2
        error('All valid depth values are identical.');
    end

    for k = 1:100
        groups = ones(size(x));
        groups(abs(x - c2) < abs(x - c1)) = 2;
        if ~any(groups == 1) || ~any(groups == 2)
            error('Two-surface grouping produced an empty group.');
        end
        next1 = mean(x(groups == 1));
        next2 = mean(x(groups == 2));
        if abs(next1 - c1) < 1e-9 && abs(next2 - c2) < 1e-9
            break;
        end
        c1 = next1;
        c2 = next2;
    end
    means = [mean(x(groups == 1)), mean(x(groups == 2))];
end

function value = interpPercentile(x, p)
    pos = 1 + (p / 100) * (numel(x) - 1);
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        value = x(lo);
    else
        value = x(lo) * (hi - pos) + x(hi) * (pos - lo);
    end
end
