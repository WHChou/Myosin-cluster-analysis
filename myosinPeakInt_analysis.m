fpath = 'E:/D/Wen-hung/03 Microscopy/Cell/20220518_Y_B_washout/';
cd(fpath);
load('myoInt.mat')

%% Read results data
sumIntYB = resultsYB.sumInt;
sumIntYBout = resultsYBout.sumInt;
sumIntYoutB = resultsYoutB.sumInt;
sumIntYoutBout = resultsYoutBout.sumInt;
totalIntYB= resultsYB.totalInt;
totalIntYBout= resultsYBout.totalInt;
totalIntYoutB= resultsYoutB.totalInt;
totalIntYoutBout= resultsYoutBout.totalInt;
modeYB = resultsYB.mode;
densityYB = resultsYB.density;
modeYBout = resultsYBout.mode;
densityYBout = resultsYBout.density;
modeYoutB = resultsYoutB.mode;
densityYoutB = resultsYoutB.density;
modeYoutBout = resultsYoutBout.mode;
densityYoutBout = resultsYoutBout.density;
sigmaYB= resultsYB.sigma;
sigmaYBout= resultsYBout.sigma;
sigmaYoutB= resultsYoutB.sigma;
sigmaYoutBout= resultsYoutBout.sigma;

%% Plot boxplots to compare mode/density in different drug perturbations
plotBoxPlot(modeYB, modeYBout, modeYoutB, modeYoutBout,  {'Y+B+', 'Y+B-', 'Y-B+', 'Y-B-'}, 'Mode of myosin cluster intensity (a.u.)');
plotBoxPlot(densityYB, densityYBout, densityYoutB, densityYoutBout, {'Y+B+', 'Y+B-', 'Y-B+', 'Y-B-'}, 'Myosin cluster density (1/um^2)');
plotBoxPlot(sigmaYB, sigmaYBout, sigmaYoutB, sigmaYoutBout, {'Y+B+', 'Y+B-', 'Y-B+', 'Y-B-'}, 'Sigma of myosin intensity distribution');

% Output results to .csv file 
writematrix([modeYB', modeYBout', modeYoutB', [modeYoutBout nan]'], 'mode_YB_washout.csv');
writematrix([densityYB', densityYBout', densityYoutB', [densityYoutBout nan]'], 'density_YB_washout.csv');
writematrix([sigmaYB', sigmaYBout', sigmaYoutB', [sigmaYoutBout nan]'], 'sigma_YB_washout.csv');

%% Histogram comparison
points1 = resultsYB.pointData.cell14; points2 = resultsYoutB.pointData.cell18;
barint = resultsYB.param.barint;
[phat1, ~] = mle(points1(:,6), 'distribution', 'lognormal', 'TruncationBounds', [barint, Inf]);
[phat2, ~] = mle(points2(:,6), 'distribution', 'lognormal', 'TruncationBounds', [barint, Inf]);
figure; h1 = histogram(points1(:,6), 'Normalization', 'pdf'); h1.BinWidth = 500; h1.FaceAlpha = 0.5; hold on;
h2 =  histogram(points2(:,6), 'Normalization', 'pdf'); h2.BinWidth = 500; h2.FaceAlpha = 0.5; 
xt1 = 1:h1.BinLimits(2); yt1 = lognpdf(xt1, phat1(1), phat1(2)); plot(xt1, yt1, 'b', 'LineWidth', 4); 
xt2 = 1:h2.BinLimits(2); yt2 = lognpdf(xt2, phat2(1), phat2(2)); plot(xt2, yt2, 'Color', [0.828, 0.258, 0.0313], 'LineWidth', 4); 
xlim([0 20000]); xlabel('Myosin intensity (a.u.)'); ylabel('PDF'); set(gca, 'FontSize', 20)
legend({'Y+B+', 'Y-B+', '', ''})

%% Scaling of large myosin stack sizes
barint = 400;
numCell = 20;
maxRange = 25000;
binSize = 1000;
edges = 0:binSize:maxRange;
binPos = binSize/2:binSize:(maxRange-binSize/2);
binVal = zeros(numCell, length(binPos));
for i = 1:numCell
    valName = ['cell' int2str(i)];
    points1 = resultsBB.pointData.(valName);
    h = histogram(points1(:,6), edges);
    binVal(i, :) = h.Values;
    binEdge = h.BinEdges;
    close gcf;
end
id = find(binPos>10000&binPos<20000);
y = mean(binVal, 1);
Y = log(y(id)');
x = binPos;
X = [ones(length(x(id)), 1), log(x(id))'];
[b, ~, r, ~, stats] = regress(Y, X);
figure; loglog(binPos, mean(binVal, 1), 'b.', 'MarkerSize', 20); hold on
errorbar(binPos, mean(binVal, 1), std(binVal, 1), 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1)
plot(exp(X), exp(b(1)+b(2).*X), 'r', 'LineWidth', 2)
xlim([barint, maxRange]); xlabel('Intensity (a.u.)'); ylabel('Count')
set(gca, 'FontSize', 20); set(gca, 'LineWidth', 1)

%% Plot myosin intensity over time (for live-cell imaging)
dt = 30;
sumInt = results.sumInt;
totalInt = results.totalInt;
modeLN = results.mode;
punctaDensity = results.density;
% Sum pixel intensity
figure; yyaxis left; plot((1:length(results.sumInt))*dt, results.sumInt, '.', 'MarkerSize', 20)
xlabel('Time (s)'); ylabel('Total myosin pixel intensity (a.u.)'); 
%set(gca, 'FontSize', 20); box off
% Sum peak intensity
%figure; 
yyaxis right; plot((1:length(results.sumInt))*dt, results.totalInt, '.', 'MarkerSize', 20)
xlabel('Time (s)'); ylabel('Total myosin cluster intensity (a.u.)'); 
set(gca, 'FontSize', 20); box off

% Mode of intensity
t0 = 11;
figure; yyaxis left; plot(((1-t0):(length(results.sumInt)-t0))*dt, results.mode, '.', 'Color', [0, 0.4470, 0.7410], 'MarkerSize', 20); hold on;
plot(((1-t0):(length(results.sumInt)-t0))*dt, movmean(results.mode, 3), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Mode of intensity distribution (a.u.)'); 
%set(gca, 'FontSize', 20); box off
% Density of puncta
%figure; results.density(17)=nan;
yyaxis right; plot(((1-t0):(length(results.sumInt)-t0))*dt, results.density, '.', 'Color', [0.8500, 0.3250, 0.0980], 'MarkerSize', 20); hold on;
plot(((1-t0):(length(results.sumInt)-t0))*dt, movmean(results.density, 3, 'omitnan'), '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Density of myosin clusters (1/um^2)'); 
set(gca, 'FontSize', 20); box off
% Sigma (tail of distribution) results.sigma(17)=nan;
figure; plot(((1-t0):(length(results.sumInt)-t0))*dt, results.sigma, '.', 'Color', [0.9290, 0.6940, 0.1250], 'MarkerSize', 20); hold on;
plot(((1-t0):(length(results.sumInt)-t0))*dt, movmean(results.sigma, 3, 'omitnan'), '-', 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Sigma of intensity distribution'); 
set(gca, 'FontSize', 20); box off
