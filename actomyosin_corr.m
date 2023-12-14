%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colocalization analysis of myosin cluster size and F-actin intensity
% 
% This script consists of two parts.
%
% Part 1: Correlation of every myosin cluster with the underlying F-actin intensity.
% Prerequisite: user has run myosin cluster intensity analysis and saved as
% .mat file. 
%
% Part 2: Correlation of myosin clusters on select F-actin bundles
% Prerequisite: One must first mark ROIs in ImageJ and export coordinate 
% information using getSelectionCoordinates.ijm. Requires isOnfiber.m and 
% linearFit.m.
%
% Note: This script is meant to be executed line by line. Each segment
% should be re-run for different conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Correlation for every myosin cluster  %%
% Set up parameters
fpath = '/Volumes/GARDEL LAB/D/Wen-hung/03 Microscopy/Cell/20211025/';  % Folder that contains image and analysis files
pxsz = 6.45/60; % pixel size in microns
cd(fpath);

numCell = 18;
basename = '20221003_IF_siNT_U2OS_NMIIA_mScar_pgipz_x_actin_cell'; 
load('myoInt_refined2.mat');    % Load analysis results of myosin cluster sizes
m_array_NT = zeros(1, numCell);
b_array_NT = zeros(1, numCell);
r_sq_array_NT = zeros(1, numCell);
colocMat_NT = [];

% Run each condition separately with different variable names
for ii = 1:numCell
    % Read image files
    formatSpec = '%02d';
    imActin = tiff_stack_read([fpath basename num2str(ii, formatSpec) '_w1642-zyla.TIF']);
    imMyosin = tiff_stack_read([fpath basename num2str(ii, formatSpec) '_w2561-zyla.TIF']);

    % Image processing
    %Find focal plane based on actin image
    zInt = reshape(sum(sum(imActin)), 1, size(imActin, 3));
    [~, focId] = max(zInt);
    %Top-hat background subtraction 
    im_foc = stackZproj(imMyosin,focId);
    im_foc_ac = imActin(:,:,focId);
    se = strel('disk', 50);
    im_bgsub = imtophat(im_foc, se);
    im_bgsub_ac = imtophat(im_foc_ac, se);

    % Colocalization of myosin peak intensity
    varName = ['cell' num2str(ii)];
    points = resultsNT.pointData.(varName);
    colocMat2 = zeros(length(points), 3);
    for k = 1:length(points)
       colocMat2(k, 1) = im_bgsub_ac(round(points(k,2)), round(points(k,1)));
       colocMat2(k, 2) = points(k,6);
       colocMat2(k, 3) = ii;
    end
    colocMat_NT = [colocMat_NT; colocMat2];
    y = colocMat2(:,2);
    x0 = ones(length(y), 1);
    x1 = colocMat2(:,1);
    X = [x0, x1];
    [b,bint,r,rint,stats] = regress(y,X);
    r_sq_array_NT(ii) = stats(1);
    m_array_NT(ii) = b(2);
    b_array_NT(ii) = b(1);
    %%% Uncomment this section to show correlation plots for each cell %%%
    %color = [0 0.5 1];
    %figure; 
    %plot(colocMat2(:,1), colocMat2(:,2), '.', 'MarkerSize', 10, 'Color', color); hold on;
    %plot(x1, b(1)+b(2)*x1, 'm', 'LineWidth', 2)
    %xlabel('Actin pixel intensity (a.u.)'); ylabel('Myosin stack intensity (a.u.)')
    %set(gca, 'FontSize', 20); xlim([0 5000])
    %pause(1); close;
end
plotBoxPlot(r_sq_array_NT, r_sq_array_AC, {'siNT', 'siACTN1'}, 'R^2 of actomyosin correlation')
plotBoxPlot(m_array_NT, m_array_AC, {'siNT', 'siACTN1'}, 'Slope of actomyosin correlation')

%% Myosin intensity on select F-actin bundles %%
% Load files
fpath = '/Volumes/GARDEL LAB/D/Wen-hung/03 Microscopy/Cell/20211025/';  % Folder that contains image and analysis files
cd(fpath)

basename = '20211025_IF_NMIIA_mScar_pMLC_x_actin_50nM_LatA_cell';
formatSpec = '%02d';
load('myoInt_refined_4.mat');  % load myosin cluster analysis results
allIntWT1025 = resultsWT.allInt;
minIntWT1025 = min(allIntWT1025(:,2));
allIntLA1025 = resultsLA.allInt;
minIntLA1025 = min(allIntLA1025(:,2));
allIntWT0304 = resultsWT.allInt;
minIntWT0304 = min(allIntWT0304(:,2));
allIntBB0304 = resultsBB.allInt;
minIntBB0304 = min(allIntBB0304(:,2));
allIntNT = resultsNT.allInt;
minIntNT = min(allIntNT(:,2));
allIntAC = resultsAC.allInt;
minIntAC = min(allIntAC(:,2));

% Set up variables to store results
r_sq_WT0304 = [0 0 0 0 0];
m_WT0304 = [0 0 0 0 0];
allOnFiberWT0304 = [];
allDensityWT0304 = [];

% Loop through all cells which have F-actin bundles marked
count = 1;
for ii=[1,4,5,7,9]
    % Read files
    imActin = tiff_stack_read([fpath basename num2str(ii, formatSpec) '_w1642-zyla.TIF']);
    imMyosin = tiff_stack_read([fpath basename num2str(ii, formatSpec) '_w2561-zyla.TIF']);
    fibers = readmatrix(['SF_LatA_cell' num2str(ii, formatSpec) '.csv']);  % csv file that has F-actin bundle coordinates
    sf_id = unique(fibers(:,3));
    varName = ['cell' int2str(ii)];
    points1 = resultsLA.pointData.(varName);
    
    %Find focal plane based on actin image
    zInt = reshape(sum(sum(imActin)), 1, size(imActin, 3));
    [~, focId] = max(zInt);
    % Top-hat background subtraction 
    im_foc = stackZproj(imMyosin,focId);
    im_foc_ac = stackZproj(imActin, focId);
    se = strel('disk', 50);
    im_bgsub = imtophat(im_foc, se);
    im_bgsub_ac = imtophat(im_foc_ac, se);
    
    [onFiberLA1025, SFL] = isOnFiber(points1, fibers, im_bgsub_ac);
    
    %%% Uncomment this section to visualize clusters on SF segments
    %figure; imshow(imcomplement(im_foc_ac), []); hold on
    %for jj = 0:length(sf_id)    
    %   plot(fibers(fibers(:,3)==jj, 1), fibers(fibers(:,3)==jj, 2), 'Color', 'm', 'LineWidth', 3)
    %end
    %scatter(onFiberWT(:,1), onFiberWT(:,2), 10, 'g', 'filled')
    %plot([770, 863.02], [1166, 1166], 'k', 'LineWidth', 8)
    
    % Myosin density on stress fibers
    density = zeros(1, size(unique(fibers(:,3)), 1));
    for jj = 1:size(unique(fibers(:,3)), 1)   % For each stress fiber
        L = SFL(jj)*0.1075;
        numPuncta = size(onFiberWT0304(onFiberWT0304(:,8)==jj-1, :), 1);
        density(jj) = numPuncta/L;
    end
    
    % Intensity correlation
    %[b4,~,~,~,stats4] = linearFit(onFiberAC(:,6), onFiberAC(:,7));
    %[b2,~,~,~,stats2] = linearFit(onFiberLA(:,6), onFiberLA(:,7));
    %[b1,~,~,~,stats1] = linearFit(onFiberWT(:,6), onFiberWT(:,7));
    %[b3,~,~,~,stats3] = linearFit(onFiberNT(:,6), onFiberNT(:,7));
    %[b5,~,~,~,stats5] = linearFit(onFiberWT0304(:,6)/minIntWT0304, onFiberWT0304(:,7));
    %[b6,~,~,~,stats6] = linearFit(onFiberBB0304(:,6)/minIntBB0304, onFiberBB0304(:,7));
    [b7,~,~,~,stats7] = linearFit(onFiberWT1025(:,6)/minIntWT1025, onFiberWT1025(:,7));
    [b8,~,~,~,stats8] = linearFit(onFiberLA1025(:,6)/minIntLA1025, onFiberLA1025(:,7));
    %[b9,~,~,~,stats9] = linearFit(onFiberNT(:,6)/minIntNT, onFiberNT(:,7));
    %[b10,~,~,~,stats10] = linearFit(onFiberAC(:,6)/minIntAC, onFiberAC(:,7));
    r_sq_WT0304(count) = stats5(1);
    m_WT0304(count) = b5(2);
    count = count+1;
    allOnFiberWT0304 = [allOnFiberWT0304; [onFiberWT0304 ones(length(onFiberWT0304), 1)*ii]];
    allDensityWT0304 = [allDensityWT0304; [density' ones(length(density), 1)*ii]];
    
    %%% Uncomment this section to plot actomyosin correlation
    %figure; scatter(onFiberWT1025(:,7), onFiberWT1025(:,6)/minIntWT1025, 40, [0 150/255 1], 'Filled');
    %xt = linspace(min(onFiberWT1025(:,7)), max(onFiberWT1025(:,7)), 20);
    %hold on; plot(xt, b7(1)+b7(2)*xt, '-', 'Color', [45/255 123/255 179/255], 'LineWidth', 3)
    %scatter(onFiberLA1025(:,7), onFiberLA1025(:,6)/minIntLA1025, 40, [1 1 0], 'Filled');
    %xt = linspace(min(onFiberLA1025(:,7)), max(onFiberLA1025(:,7)), 20);
    %hold on; plot(xt, b8(1)+b8(2)*xt, '-', 'Color', [179/255 178/255 18/255], 'LineWidth', 3)
    %xlabel('F-actin bundle intensity (a.u.)'); ylabel('Myosin cluster size (x N_{min})')
    %set(gca, 'FontSize', 20); legend({'Ctrl', '', 'LatA', ''})
end
plotBoxPlot(r_sq_WT0304, r_sq_BB0304, {'Ctrl', 'Bleb'}, 'R^2 of actomyosin correlation')
plotBoxPlot(m_WT0304, m_BB0304, {'Ctrl', 'Bleb'}, 'Slope of actomyosin correlation')
plotBoxPlot(r_sq_WT1025, r_sq_LA1025, {'Ctrl', 'LatA'}, 'R^2 of actomyosin correlation')
plotBoxPlot(m_WT1025, m_LA1025, {'Ctrl', 'LatA'}, 'Slope of actomyosin correlation')
plotBoxPlot(r_sq_NT, r_sq_AC, {'siNT', 'siACTN1'}, 'R^2 of actomyosin correlation')
plotBoxPlot(m_NT, m_AC, {'siNT', 'siACTN1'}, 'Slope of actomyosin correlation')