function smRm = smooth_2d_Ca_ratemap(ratemap, binSz, varargin)
% function smRm = smooth_2d_Ca_ratemap(ratemap, binSz, varargin)
% 
% PURPOSE:
%   Smooth 2-d ratemap with a gaussian kernel.
% 
% INPUT:
%   ratemap = unsmooth ratemap
%   binSz = spatial bin size for ratemap, in cm
%   gWinSz = (optional) size of gaussian smoothing window, in cm
%       - Default = 25 cm
%   gWinStd = (optional) standard deviation of gaussian smoothing window, in cm
%       - Default = gWinSz / 6
% 
% OUTPUT:
%   smRm = smoothed ratemap
% 
% MMD - based on code by JBT
% 01/2024
% Colgin Lab

if nargin == 2
    gWinSz = 25; %cm
else
    gWinSz = varargin{1};
end

if nargin < 4
    gWinStd = gWinSz / 6;
else 
    gWinStd = varargin{2};
end

%convert window size to # of bins based on binSz
gWinSz = round(gWinSz/binSz);
gWinStd = round(gWinStd/binSz);
gKrnl = fspecial('gaussian', gWinSz, gWinStd);

smRm = nanconv(ratemap, gKrnl);
smRm(isnan(ratemap)) = NaN; %keep unvisisted parts unvisited

end %function