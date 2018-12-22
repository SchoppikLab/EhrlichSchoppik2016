% a function to return binned xy data
% takes x, associate y values, and the number you want in each bin
% returns vectors of the mean x, mean associate y, and x and y standard
% deviation.

function [xout, yout, xstdout, ystdout] = makeEvenHistogram(x,y,numPerBin)

% pad x and y so that the last bin works OK
if (mod(length(x),numPerBin))
pad = nan(numPerBin - mod(length(x),numPerBin),1);

x = [x(:);pad];
y = [y(:);pad];
end

% make sure that X and Y are the same length
if numel(x) ~= numel(y)
  disp('X and Y should be the same length');
    return;
end


lastBin = floor(length(x)/numPerBin);

[~, dex] = sort(x);
dex = reshape(dex(1:(numPerBin*lastBin)),numPerBin,lastBin);


% % % xout = nanmean(x(dex));
% % % xstdout = nanstd(x(dex));
xout = mean(x(dex));
xstdout = std(x(dex));

yByCols = y(dex);

% % % yout = nanmean(yByCols); 
% % % ystdout = nanstd(yByCols);

yout = mean(yByCols);
ystdout = std(yByCols);

end

