%% Mirror-displacement calibration using SFG delay scans

clc; close all;

% Scan settings
depthSteps_um = [10, 20, 30, 40, 50, 75, 100];
mirrorPositions_um = [0, cumsum(depthSteps_um)];
Nrepeats = 50;

delayRange = 632.5:0.1:636.5;
odlSettle_s = 0.1;
stageSettle_s = 1.0;
moveTimeout_ms = 60000;
h = 0.75;
delayToDepth_um = 299.792458 / 2;

% Run selection
maxAbsStepErrorFactor = 3;
maxRunStdErrorFactor = 3;
enableRunRejection = true;

% Display settings
fontSize_axes = 20;
fontSize_labels = 24;
lineWidth = 1.8;

if ~exist('kdc', 'var') || ~exist('mdl', 'var') || ~exist('pm', 'var')
    error('Initialize the KDC101 stage (kdc), delay line (mdl), and SFG detector (pm) first.');
end

nPos = numel(mirrorPositions_um);
nSteps = numel(depthSteps_um);
nDelay = numel(delayRange);
allTraces_raw = nan(Nrepeats, nPos, nDelay);
centerDelay_ps = nan(Nrepeats, nPos);
leftCrossing_ps = nan(Nrepeats, nPos);
rightCrossing_ps = nan(Nrepeats, nPos);
measuredStep_um = nan(Nrepeats, nSteps);

kdc.MoveTo(System.Decimal(0), moveTimeout_ms);
pause(stageSettle_s);
zeroStagePos_mm = str2double(char(kdc.Position.ToString()));

for rep = 1:Nrepeats
    fprintf('Repeat %d / %d\n', rep, Nrepeats);

    for p = 1:nPos
        targetPosition_mm = zeroStagePos_mm + mirrorPositions_um(p) / 1000;
        kdc.MoveTo(System.Decimal(targetPosition_mm), moveTimeout_ms);
        pause(stageSettle_s);

        mdl.setAbsoluteDelay(delayRange(1));
        pause(5);
        trace = nan(1, nDelay);

        for k = 1:nDelay
            mdl.setAbsoluteDelay(delayRange(k));
            pause(odlSettle_s);
            fprintf(pm, 'MEAS:SCAL:FREQ?');
            trace(k) = str2double(fscanf(pm));
        end

        allTraces_raw(rep, p, :) = trace;
        [centerDelay_ps(rep, p), leftCrossing_ps(rep, p), ...
            rightCrossing_ps(rep, p)] = crossingCenter(delayRange, trace, h);

        fprintf('Mirror %.0f um, measured center %.4f ps\n', ...
            mirrorPositions_um(p), centerDelay_ps(rep, p));
    end

    measuredStep_um(rep, :) = diff(centerDelay_ps(rep, :)) * delayToDepth_um;
end

stepError_um = measuredStep_um - depthSteps_um;
runHasNaN = any(~isfinite(centerDelay_ps), 2) | any(~isfinite(measuredStep_um), 2);
runStepErrorTooLarge = any(abs(stepError_um) > ...
    maxAbsStepErrorFactor * depthSteps_um, 2);
runStdErrorTooLarge = std(stepError_um, 0, 2, 'omitnan') > ...
    maxRunStdErrorFactor * min(depthSteps_um);

if enableRunRejection
    badRuns = runHasNaN | runStepErrorTooLarge | runStdErrorTooLarge;
else
    badRuns = runHasNaN;
end
useRuns = find(~badRuns);
if isempty(useRuns)
    error('No valid calibration runs remain. Inspect the recorded traces.');
end

meanMeasuredStep_um = mean(measuredStep_um(useRuns, :), 1);
stdMeasuredStep_um = std(measuredStep_um(useRuns, :), 0, 1);
semMeasuredStep_um = stdMeasuredStep_um / sqrt(numel(useRuns));

fprintf('\nRuns included: %s\n', mat2str(useRuns.'));
fprintf('Actual (um)    Measured (um)    STD (um)\n');
for p = 1:nSteps
    fprintf('%10.1f    %13.3f    %8.3f\n', depthSteps_um(p), ...
        meanMeasuredStep_um(p), stdMeasuredStep_um(p));
end

figure('Color', 'w', 'Name', 'Mirror displacement calibration');
hold on;
errorbar(depthSteps_um, meanMeasuredStep_um, stdMeasuredStep_um, 'o', ...
    'LineWidth', lineWidth, 'MarkerSize', 8, 'CapSize', 10);
axisMax = 1.1 * max([depthSteps_um, meanMeasuredStep_um + stdMeasuredStep_um]);
plot([0 axisMax], [0 axisMax], '--k', 'LineWidth', 1.5);
xlabel('Actual mirror displacement ($\mu$m)', ...
    'Interpreter', 'latex', 'FontSize', fontSize_labels);
ylabel('Measured depth shift ($\mu$m)', ...
    'Interpreter', 'latex', 'FontSize', fontSize_labels);
legend('Measured depth', 'Ideal response', 'Location', 'northwest');
axis equal;
xlim([0 axisMax]); ylim([0 axisMax]);
grid on; box on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontSize_axes);

repToPlot = useRuns(1);
referenceDelay_ps = centerDelay_ps(repToPlot, 1);
figure('Color', 'w', 'Name', 'Measured SFG delay traces');
hold on;
for p = 1:nPos
    trace = squeeze(allTraces_raw(repToPlot, p, :)).';
    if all(isfinite(trace)) && max(trace) > 0
        plot((delayRange - referenceDelay_ps) * delayToDepth_um, ...
            trace / max(trace), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%.0f $\mu$m', mirrorPositions_um(p)));
    end
end
yline(h, '--k', sprintf('h = %.2f', h), 'HandleVisibility', 'off');
xlabel('Relative depth ($\mu$m)', ...
    'Interpreter', 'latex', 'FontSize', fontSize_labels);
ylabel('Normalized SFG signal', ...
    'Interpreter', 'latex', 'FontSize', fontSize_labels);
legend('Interpreter', 'latex', 'Location', 'best');
grid on; box on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontSize_axes);

fileName = ['Mirror_Steps_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.mat'];
save(fileName, 'allTraces_raw', 'delayRange', 'depthSteps_um', ...
    'mirrorPositions_um', 'centerDelay_ps', 'leftCrossing_ps', ...
    'rightCrossing_ps', 'measuredStep_um', 'stepError_um', 'useRuns', ...
    'badRuns', 'h', 'Nrepeats', 'zeroStagePos_mm');
fprintf('Saved %s\n', fileName);

function [center, left, right] = crossingCenter(t, trace, h)
    center = NaN; left = NaN; right = NaN;
    if any(~isfinite(trace)) || max(trace) <= 0
        return;
    end
    trace = trace(:).' / max(trace);
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
