function [wDist,xaxis] = WeightedProportion(x,xmin,xmax,sigma)
    xaxis = xmin:sigma/5:xmax;
    delta = repmat(x,1,length(xaxis))-repmat(xaxis,length(x),1);
    W = exp(-0.5*delta.*delta/sigma^2);
    W = W./repmat(sum(W,2),1,size(W,2)); %normalize each kernel so the ones on the edge have same contribution
    wDist = sum(W)./sum(sum(W));


end %function