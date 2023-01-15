%% Specify file path and movie name
fpath = 'E:\D\Wen-hung\03 Microscopy\Cell\20220301\';
fname = '20220301_LCI_50uM_blebbistatin_washout_NMIIA_mScar_cell3.tif';
cd(fpath)

%% Background subtraction
% Read file
myosinStack = tiff_stack_read([fpath fname]);
%%% Uncomment if movie is .czi file
%stack = czi_stack_readCZT([fpath fname]);
%myosinStack = squeeze(stack.c1);

% Estimate background
cellMask = struct();
bkgInt = zeros(1, size(myosinStack, 3));
for ii = 1:size(myosinStack, 3)
    imMyosin = myosinStack(:,:,ii);

    if mod(ii, 5) == 1
        mask = createCellMask(imMyosin, 'manual');  % Use 'auto' for Patrick's code
        varName = ['t' int2str(ii)];
        cellMask.(varName) = mask;
    end
    %cellArea = sum(sum(mask))*pxsz^2;
    bkgMask = ~mask;   % Isolate backgroubnd pixels
    notcell = imMyosin .* uint16(bkgMask);
    bkgInt(ii) = sum(sum(notcell))/sum(sum(bkgMask));
end
meanBkg = mean(bkgInt);
myosinStackSub = myosinStack - meanBkg;
myosinStackSub(myosinStackSub<0) = 0;
%tiff_stack_write(myosinStackSub, '561_align_crop_bgsub.tif');

%% Photobleaching correction: total myosin intensity should remain constant
% Establish photobleaching curve
myosinInt = zeros(1, size(myosinStack, 3));
for ii = 1:size(myosinStack, 3)
    if mod(ii,5)==1
        varName = ['t' int2str(ii)];
        mask = cellMask.(varName);
    end
    myosinInt(ii) = sum(sum(myosinStackSub(:,:,ii) .* uint16(mask)));%sum(sum(myosinStack(:,:,ii)));
end
intensity = myosinInt ./ myosinInt(1);
t = 1:size(myosinStack, 3);

% Fit two-exponential model Ae^(-bt)+(1-A)e^(-ct)
fun = @(par,t)par(1)*(exp(-par(2)*t))+(1-par(1))*(exp(-par(3)*t));
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
[par, resnorm, res] = lsqcurvefit(fun,[0 0 0], t, intensity, [], [], options);
figure; plot(t, intensity, 'r.', t, fun(par,t), 'b', 'MarkerSize', 10, 'LineWidth', 2)
pause(3); close;

% Correct using established photobleaching curve 
myosinStackCorr = zeros(size(myosinStack, 1), size(myosinStack, 2), size(myosinStack, 3));
for ii = 1:size(myosinStack, 3)
    imMyosin = myosinStackSub(:,:,ii); 
    myosinStackCorr(:,:,ii) = imMyosin./fun(par, ii);
    imwrite(uint16(myosinStackCorr(:,:,ii)), '20220301_LCI_50uM_blebbistatin_washout_NMIIA_mScar_cell3_corr.tif', 'WriteMode', 'append')
end