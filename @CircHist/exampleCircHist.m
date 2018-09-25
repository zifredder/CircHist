%% Using the CircHist Class
%% Plot Distribution Data
% Generate a noisy sample (von Mises distribution with |theta| == 90 deg).
rng default
sDist = mod(rad2deg(circ_vmrnd(pi/2, 2, 100)), 360); % generate sample, convert to deg
nBins = 36; % number of bins, makes bin size of 10 deg
%%
% Plot the circular histogram:
obj1 = CircHist(sDist, nBins);
%%
% Adjust appearance:
obj1.colorBar = 'k'; 
obj1.avgAngH.LineStyle = '--';
obj1.avgAngH.LineWidth = 1;
obj1.colorAvgAng = [.5 .5 .5];
% remove offset between bars and plot-center
rl = rlim;
obj1.setRLim([0, rl(2)]);
% draw circle at r == 0.5 (where r == 1 would be the outer plot edge)
rl = rlim;
obj1.drawCirc((rl(2) - rl(1)) /2, '--b', 'LineWidth', 2)
obj1.scaleBarSide = 'right';
obj1.polarAxs.ThetaZeroLocation = 'right';
obj1.setThetaLabel('Direction', 'bottomleft');
% draw resultant vector r as arrow
delete(obj1.rH)
obj1.drawArrow(obj1.avgAng, obj1.r * range(rl), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'r')
% Change theta- and rho-axis ticks
obj1.polarAxs.ThetaAxis.MinorTickValues = [];
thetaticks(0:90:360);
rticks(0:4:20);
obj1.drawScale; % update scale bar
%% Plot Multi-Sample Distribution
% Generate another noisy sample with a different distribution-width |kappa|.
rng default
s2Dist = mod(rad2deg(circ_vmrnd(pi/2, 1.5, 100)), 360);
sCell = {sDist, s2Dist}; % pack both samples into a cell-array
figure
CircHist(sCell, nBins);
%% Combine Multiple Histograms in One Figure
% Create subplot, note that the created |axes| must be a |polaraxes|.
fH = figure;
subAx1 = subplot(1, 2, 1, polaraxes);
subAx2 = subplot(1, 2, 2, polaraxes);
obj2 = CircHist(sDist, nBins, 'ax', subAx1);
obj3 = CircHist(s2Dist, nBins, 'ax', subAx2);
% Make rho-axes equal for both diagrams
maxRho = max([max(rlim(subAx1)), max(rlim(subAx2))]);
newLimits = [min(rlim(subAx1)), maxRho];
obj2.setRLim(newLimits);
obj3.setRLim(newLimits);
% Adjust figure-window size
drawnow
fH.Position(3) = 1100; % width
fH.Position(4) = 500; % height
%% Plot Already-Binned Data
% Bin the generated multi-sample distribution before plotting.
edges = 0:10:360;
histData = histcounts(mod([sDist; s2Dist], 360), edges);
% Note that |edges| can be omitted because the number of bins results from the number of
% data points in |histData|, but that |'dataType'| must be specified as |'histogram'|.
figure
CircHist(histData, 'dataType', 'histogram');
%% Axial Data
% Copy the von Mises data with an offset of 180 deg and a little bit of noise to generate
% an axial, bimodal distribution.
rng default
noise = (rand(size(sDist)) - 0.5) * 10;
sAxial = [sDist; sDist + 180 + noise];
%%
% Call |CircHist| with |'areAxialData'| specified as |true|.
figure
CircHist(sAxial, nBins, 'areAxialData', true);
%%
% Note that now the average angle is indicated by an axis that halves the diagram at this
% angle.
%% Draw Arrows
rng default
arrowLen = randn(numel(sDist), 1); % random arrow lengths
arrowLen = arrowLen / max(arrowLen);
arrowLen = arrowLen + abs(min(arrowLen));
figure
obj4 = CircHist([1, 2], 36); % dummy data
delete([obj4.avgAngH; obj4.avgAngCiH(:); obj4.barH(:); obj4.rH]); % remove dummy data
title('');
obj4.scaleBar.Label.String = 'Vector length';
obj4.polarAxs.ThetaAxis.MinorTickValues = [];
thetaticks(0:90:360);
arrowH = obj4.drawArrow(sDist, arrowLen);
%%
% Change visual properties and add another arrow.
set(arrowH, 'HeadStyle', 'plain', 'HeadWidth', 3)
% Draw a single arrow that ends at the outer plot edge
avgAng = circ_mean(deg2rad(sDist), arrowLen); % average angle, weighted by arrow length
obj4.drawArrow(rad2deg(avgAng), [], 'Color', 'r', 'LineWidth', 3)
%%
drawnow % (else, the last figure is not shown in the published version for some reason)
%% Enable Tab Auto-Completion for Object Construction
% If |functionSignatures.json| is located in the same directory as the |@CircHist| folder,
% Name-Value pairs of the object-constructor call can be auto-completed as it is the case
% for builtin MATLAB functions. See also:
% <https://mathworks.com/help/matlab/matlab_prog/customize-code-suggestions-and-completions.html>
%
% <<tab-auto.png>>
%