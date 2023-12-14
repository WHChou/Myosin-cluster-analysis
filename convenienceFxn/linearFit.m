function [b,bint,r,rint,stats] = linearFit(ydata, xdata)
    x0 = ones(length(ydata), 1);
    X = [x0, xdata];
    [b,bint,r,rint,stats] = regress(ydata,X);
end