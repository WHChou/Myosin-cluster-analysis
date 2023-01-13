%% Parameters to load data
% Specify file path
fpath = 'E:\D\Wen-hung\03 Microscopy\Cell\20221014\';  
cd(fpath);

% Specify the number of images & image name
numCell = 16;
basename = '20211020_IF_NMIIA_mScar_pMLC_x_actin_ctrl_cell'; 
pxsz = 6.45/60; 

%% Set up arrays to store results
pointData = struct();        
cellMask = struct();
cellArea = zeros(1, numCell);
sumInt = zeros(1, numCell);
totalInt = zeros(1, numCell);
numPuncta = zeros(1, numCell);
density = zeros(1, numCell);
modeLN = zeros(1, numCell);
med = zeros(1, numCell);
sigma = zeros(1, numCell);

%% Start analysis
for ii = 1:numCell
    ii
    % Read file
    formatSpec = '%02d';
    imActin = tiff_stack_read([fpath basename num2str(ii, formatSpec) '_w1642-zyla.TIF']);
    imMyosin = tiff_stack_read([fpath basename num2str(ii, formatSpec) '_w2561-zyla.TIF']);

    %%% Pre-processing
    % Find focal plane based on actin image & perform projection
    zInt = reshape(sum(sum(imActin)), 1, size(imActin, 3));
    [~, focId] = max(zInt);
    im_foc = stackZproj(imMyosin,focId);
    im_foc_ac = stackZproj(imActin,focId);

    % Top-hat background subtraction
    se = strel('disk', 50); %offsetstrel('ball', 50, 3);
    im_bgsub = imtophat(im_foc, se);

    % Estimate background variance (2*variance is the signal floor)
    bkgMask = ~createCellMask(imActin(:,:,focId), 'auto1');
    crop_coorx = [size(im_foc, 1)/4, size(im_foc, 1)*3/4];  % Assume center of frame has uniform bkg
    crop_coory = [size(im_foc, 2)/4, size(im_foc, 2)*3/4];  % So crop central region and calculate bkg
    crop_m = im_bgsub(crop_coorx(1):crop_coorx(2), crop_coory(1):crop_coory(2));
    crop_mask = bkgMask(crop_coorx(1):crop_coorx(2), crop_coory(1):crop_coory(2));
    im_bg = crop_m .* uint16(crop_mask);
    bkgIntArr = im_bg(im_bg~=0);
    bkgVar = std(double(bkgIntArr(bkgIntArr<100)));

    % Find myosin puncta
    rad = 3;
    barint = 100;   % min feature intensity
    barrg = 4.75;   % radius of gyration: filter out features that are too "spread out"
    barcc = 0.48;   % eccentricity
    IdivRg = 180;   % intensity/radius of gyration: filter out features whose intensity aren't distributed evenly
    masscut = 0;
    Imin = ceil(2*bkgVar);   % min pixel intensity to be considered a puncta
    [~, r, rej] = mpretrack_wc(im_bgsub_denoise, rad, barint, barrg, barcc, IdivRg, masscut, Imin);
    %%% Uncoment to check if feature finding parameters are reasonable
    %figure; imshow(im_bgsub, []); hold on;
    %plot(r(:,1), r(:,2), 'r.', 'MarkerSize', 10)
    %plot(rej(:,1), rej(:,2), 'g.', 'MarkerSize', 10)

    % Create cell mask & filter out features outside of the cell
    mask = createCellMask(im_bgsub, 'auto1');  % Use 'auto' for Patrick's code

    points = [];
    for jj = 1:length(r)
        x = round(r(jj, 2));
        y = round(r(jj, 1));
        if mask(x, y)~=0 %& mask_nuc(x, y)==0
            points = [points; r(jj, :)];
        end
    end
    
   %% Filter out points that are too close, save figure
    [points1, Erej] = filtClosePoints(points, 2);
    
%%% Uncomment to overlay identified clusters on original image
%     plotPuncta(im_bgsub_denoise, points1, rej, Erej, 8);       
%     filename = ['puncta map/20221014_LCI_YoutBleb_cell3_puncta_t' int2str(ii) '.png'];
%     f = gcf;
%     exportgraphics(f, filename, 'Resolution', 600);
%     close;

%%% Uncomment to overlay color-coded clusters on original image
%     plotColorPuncta(im_bgsub, points1, 2000, 1, 1)
%     filename = ['20211020_ctrl_peaks_color_cell' int2str(ii) '.png'];
%     f = gcf;
%     exportgraphics(f, filename, 'Resolution', 300);
%     close;
    
    %% Process intensity distribution
    [phat, ~] = mle(points1(:,6), 'distribution', 'lognormal', 'TruncationBounds', [0, Inf]); %barint
%%% Uncomment to plot intensity histogram and lognormal fit
%     figure; h = histogram(points1(:,6), 'Normalization', 'pdf'); h.BinWidth = 500;
%     xt = 1:h.BinLimits(2);
%     yt = lognpdf(xt, phat(1), phat(2));
%     % % yt = pdf(pd, xt);
%     hold on; plot(xt, yt, 'r', 'LineWidth', 2); xlim([0 10000])
%     xlabel('Intensity (a.u.)'); ylabel('PDF'); set(gca, 'FontSize', 20)
%     filename = ['histogram/20221014_LCI_YoutBleb_cell3_histogram_t' int2str(ii) '.png'];%num2str(ii, formatSpec)
%     f = gcf;
%     exportgraphics(f, filename, 'Resolution', 300);
%     close;

   %% Save results
    intensity = points1(:, 6);
    varName = ['cell' int2str(ii)];
    pointData.(varName) = points1;        
    cellMask.(varName) = mask;
    cellArea(ii) = sum(sum(mask))*pxsz^2;
    sumInt(ii) = sum(sum(im_bgsub .* uint16(mask)));
    totalInt(ii) = sum(intensity);
    density(ii) = length(points1)/cellArea(ii);
    numPuncta(ii) = length(points1);
    modeLN(ii) = exp(phat(1)-phat(2)^2);
    sigma(ii) = phat(2);    
    med(ii) = exp(phat(1));
end

%% Save results under a single struct
results = struct();
param = struct();
param.rad = rad;
param.barint = barint;
param.barrg = barrg;
param.barcc = barcc;
param.IdivRg = IdivRg;
param.masscut = masscut;
param.Imin = Imin;
results.pointData = pointData;
results.cellMask = cellMask;
results.cellArea = cellArea;
results.sumInt = sumInt;
results.totalInt = totalInt;
results.density = density;
results.numPuncta = numPuncta;
results.mode = modeLN;
results.sigma = sigma;
results.med = med;
results.param = param;
save('myoInt.mat', 'results')
